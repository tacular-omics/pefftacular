r"""Applying PEFF sequence variants to a sequence.

``\VariantSimple`` and ``\VariantComplex`` describe sequence variations the entry's
sequence line does **not** contain (spec §3.3.8, §3.3.9). This module turns a chosen
set of them into a real sequence, remaps every other annotation onto it, and enumerates
the forms an entry describes.

The arithmetic lives here; :class:`~pefftacular.VariantForm` is plain data in
``_models.py`` and :class:`~pefftacular.SequenceEntry` only delegates, as it does for
``to_proforma``.
"""

from __future__ import annotations

import dataclasses
import itertools
import math
import warnings
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Literal, overload

from pefftacular._models import (
    DisulfideBond,
    ModRes,
    ModResPsi,
    ModResUnimod,
    Processed,
    Proteoform,
    SequenceEntry,
    SequenceRange,
    VariantComplex,
    VariantForm,
    VariantSimple,
)
from pefftacular.errors import PeffError, PeffWarning

Variant = VariantSimple | VariantComplex

#: Default ceiling on the number of forms :func:`entry_variant_forms` will enumerate.
DEFAULT_FORM_LIMIT = 10_000

STOP = "*"


def _error(message: str, hint: str) -> PeffError:
    err = PeffError(message)
    err.add_note(f"hint: {hint}")
    return err


def _int_position(value: int | str, *, what: str, name: str) -> int:
    """Return *value* as an ``int``, or raise: a variant must be locatable."""
    if isinstance(value, int):
        return value
    try:
        return int(value)
    except (TypeError, ValueError):
        raise _error(
            f"{name}: {what} position {value!r} is not a number",
            "a variant can only be applied at a numeric position; '?' positions are allowed "
            "for modifications (spec section 3.3.10), not for variants",
        ) from None


def _source_span(variant: Variant, length: int, name: str) -> tuple[int, int]:
    """The 1-based inclusive range of residues *variant* replaces.

    Spec §3.3.8 requires a ``VariantSimple`` position in ``1..length``; §3.3.9 gives
    ``VariantComplex`` a start and an end that are both 1-based and inclusive.
    """
    if isinstance(variant, VariantSimple):
        pos = _int_position(variant.position, what="VariantSimple", name=name)
        if not 1 <= pos <= length:
            raise _error(
                f"{name}: VariantSimple position {pos} is outside 1..{length}",
                "spec section 3.3.8: the position must be greater than 0 and at most the protein length",
            )
        return pos, pos

    start = _int_position(variant.start_pos, what="VariantComplex start", name=name)
    end = _int_position(variant.end_pos, what="VariantComplex end", name=name)
    if not 1 <= start <= end <= length:
        raise _error(
            f"{name}: VariantComplex range {start}-{end} is outside 1..{length}",
            "spec section 3.3.9: start and end are 1-based, inclusive, and start must not exceed end",
        )
    return start, end


def _check_new_amino_acid(variant: VariantSimple, pos: int, name: str) -> None:
    """A ``VariantSimple`` substitutes one residue, or ``*`` for a stop (spec §3.3.8)."""
    new = variant.new_amino_acid
    if new == STOP:
        return
    if len(new) != 1 or not ("A" <= new <= "Z"):
        raise _error(
            f"{name}: VariantSimple at {pos}: new amino acid {new!r} is not a single letter A-Z or '*'",
            "spec section 3.3.8: a valid amino acid code (ambiguity codes included) or '*' for a stop codon",
        )


def _check_new_sequence(variant: VariantComplex, start: int, end: int, name: str) -> None:
    """A ``VariantComplex`` replaces residues with residue letters, or with nothing.

    Spec §3.3.9 gives ``(100|100|)`` for a deletion and says no character may stand for
    one, allows no regular expressions, and never gives the asterisk a meaning here: the
    nonsense codon belongs to ``VariantSimple`` (§3.3.8). So the replacement is either
    empty or amino acid letters, checked the way ``_check_new_amino_acid`` checks one.
    """
    new = variant.new_sequence
    if new == "" or all("A" <= char <= "Z" for char in new):
        return
    raise _error(
        f"{name}: VariantComplex ({start}|{end}|{new}): new sequence {new!r} is not amino acid letters A-Z",
        "spec section 3.3.9: the new sequence is residue letters, or empty for a deletion - no dash, "
        "no regular expression, and no '*' (the stop codon is a VariantSimple, section 3.3.8)",
    )


def _ordered(variants: Iterable[Variant], length: int, name: str) -> list[tuple[int, int, Variant]]:
    """Validate the variants and return them sorted by position, rejecting overlaps.

    The spec does not say how two variants over the same residue combine, so applying
    both is a conflict rather than a guess.
    """
    spans: list[tuple[int, int, Variant]] = []
    for variant in variants:
        start, end = _source_span(variant, length, name)
        if isinstance(variant, VariantSimple):
            _check_new_amino_acid(variant, start, name)
        else:
            _check_new_sequence(variant, start, end, name)
            if len(variant.new_sequence) == 1 and start == end:
                # (100|100|A) is a single-residue substitution and MUST be a VariantSimple
                # (spec section 3.3.9). Reading is permissive, so apply it and say so.
                warnings.warn(
                    f"{name}: VariantComplex ({start}|{end}|{variant.new_sequence}) is a single-residue "
                    "substitution and should be a VariantSimple (spec section 3.3.9)",
                    PeffWarning,
                    stacklevel=5,
                )
        spans.append((start, end, variant))

    spans.sort(key=lambda item: (item[0], item[1]))
    for (_start, end, variant), (next_start, _, next_variant) in itertools.pairwise(spans):
        if next_start <= end:
            raise _error(
                f"{name}: variants {_describe(variant)} and {_describe(next_variant)} both change "
                f"position {next_start}",
                "the spec does not define how overlapping variants combine; apply one per form",
            )
    for index, (start, _end, variant) in enumerate(spans[:-1]):
        if isinstance(variant, VariantSimple) and variant.new_amino_acid == STOP:
            # The protein ends here (spec section 3.3.8), so the positions the later
            # variants name do not exist in this form. Refusing keeps ``applied``
            # honest: every variant it lists is in the sequence.
            raise _error(
                f"{name}: {_describe(spans[index + 1][2])} is after the stop {_describe(variant)} at {start}",
                "a nonsense variant ends the protein (spec section 3.3.8), so nothing downstream of it "
                "can be applied in the same form; apply them in separate forms",
            )
    return spans


def _describe(variant: Variant) -> str:
    if isinstance(variant, VariantSimple):
        return f"VariantSimple({variant.position}|{variant.new_amino_acid})"
    return f"VariantComplex({variant.start_pos}|{variant.end_pos}|{variant.new_sequence})"


def _splice(
    sequence: str, spans: Sequence[tuple[int, int, Variant]]
) -> tuple[str, list[int | None], set[int], int | None]:
    """Build the variant sequence, its position map, substituted sites and any stop.

    The map is 1-based on both sides, index 0 unused, ``None`` for a residue the variant
    consumed. ``substituted`` holds the original positions whose residue a
    ``VariantSimple`` changed: a modification there is dropped (spec §3.3.10).
    """
    length = len(sequence)
    pieces: list[str] = []
    position_map: list[int | None] = [None] * (length + 1)
    substituted: set[int] = set()
    truncated_at: int | None = None
    cursor = 1  # next original position not yet emitted
    out_len = 0

    for start, end, variant in spans:
        pieces.append(sequence[cursor - 1 : start - 1])
        for pos in range(cursor, start):
            out_len += 1
            position_map[pos] = out_len

        if isinstance(variant, VariantSimple):
            if variant.new_amino_acid == STOP:
                # A nonsense mutation: the protein ends before this position and a peptide
                # must not span it (spec section 3.3.8). Nothing after it exists.
                truncated_at = start
                return "".join(pieces), position_map, substituted, truncated_at
            pieces.append(variant.new_amino_acid)
            out_len += 1
            position_map[start] = out_len
            substituted.add(start)
        else:
            pieces.append(variant.new_sequence)
            out_len += len(variant.new_sequence)
            # Every replaced residue is gone, even when the replacement is longer: the
            # new residues are not the old ones (spec section 3.3.9).
            for pos in range(start, end + 1):
                position_map[pos] = None
        cursor = end + 1

    pieces.append(sequence[cursor - 1 :])
    for pos in range(cursor, length + 1):
        out_len += 1
        position_map[pos] = out_len
    return "".join(pieces), position_map, substituted, truncated_at


def _warn_dropped(name: str, what: str, why: str) -> None:
    # Six frames out is the caller of ``SequenceEntry.apply_variants``: this helper, the
    # per-annotation remapper, ``_apply``, ``entry_apply_variants`` and the method.
    warnings.warn(f"{name}: {what} dropped: {why}", PeffWarning, stacklevel=6)


def _map_one(position: int | str, position_map: Sequence[int | None]) -> int | str | None:
    """Map one annotation position; ``None`` means the residue is gone.

    A position that is not a number (``?``, spec §3.3.10) is kept as written: it does not
    name a residue, so nothing can move it.
    """
    if not isinstance(position, int):
        try:
            numeric = int(position)
        except (TypeError, ValueError):
            return position
    else:
        numeric = position
    if not 1 <= numeric < len(position_map):
        return None
    return position_map[numeric]


def _remap_mods[Mod: (ModResUnimod, ModResPsi, ModRes)](
    mods: tuple[Mod, ...],
    position_map: Sequence[int | None],
    substituted: set[int],
    name: str,
    dropped_ids: set[int],
) -> tuple[Mod, ...]:
    """Remap modification sites, dropping those the variant removed or changed."""
    kept: list[Mod] = []
    for mod in mods:
        positions: list[int | str] = []
        for position in mod.positions:
            if isinstance(position, int) and position in substituted:
                # Spec section 3.3.10: a modified variant needs its own entry, so a
                # modification on a substituted residue cannot be carried over.
                _warn_dropped(
                    name, f"{type(mod).__name__} {mod.accession} at {position}", "the residue was substituted"
                )
                continue
            mapped = _map_one(position, position_map)
            if mapped is None:
                _warn_dropped(name, f"{type(mod).__name__} {mod.accession} at {position}", "the residue was removed")
                continue
            positions.append(mapped)
        if not positions:
            if mod.annot_id is not None:
                dropped_ids.add(mod.annot_id)
            continue
        kept.append(dataclasses.replace(mod, positions=tuple(positions)))
    return tuple(kept)


def _remap_range(
    start: int | str, end: int | str, position_map: Sequence[int | None]
) -> tuple[int | str, int | str] | None:
    """Remap both endpoints of a range, or ``None`` if either is gone.

    A range whose endpoint the variant removed is dropped rather than clipped: a clipped
    range is a claim the file does not make.
    """
    new_start = _map_one(start, position_map)
    new_end = _map_one(end, position_map)
    if new_start is None or new_end is None:
        return None
    return new_start, new_end


def _remap_processed(
    processed: tuple[Processed, ...], position_map: Sequence[int | None], name: str, dropped_ids: set[int]
) -> tuple[Processed, ...]:
    kept: list[Processed] = []
    for item in processed:
        mapped = _remap_range(item.start_pos, item.end_pos, position_map)
        if mapped is None:
            _warn_dropped(name, f"Processed {item.accession} {item.start_pos}-{item.end_pos}", "a bound was removed")
            if item.annot_id is not None:
                dropped_ids.add(item.annot_id)
            continue
        kept.append(dataclasses.replace(item, start_pos=mapped[0], end_pos=mapped[1]))
    return tuple(kept)


def _remap_pending_variants(
    entry: SequenceEntry,
    applied: Mapping[Variant, int],
    position_map: Sequence[int | None],
    name: str,
    dropped_ids: set[int],
) -> tuple[tuple[VariantSimple, ...], tuple[VariantComplex, ...]]:
    """Carry the variants that were *not* applied onto the new coordinates.

    An annotation counts as applied when it *equals* one that was applied, not when it is
    the same object: the models are frozen values, so a caller that rebuilds a variant
    instead of passing the entry's own instance means the same variant. One occurrence is
    consumed per application, so an entry that lists the same variant twice keeps the
    second as pending.
    """
    pending = Counter(applied)
    simple: list[VariantSimple] = []
    for variant in entry.variant_simple:
        if pending[variant]:
            pending[variant] -= 1
            continue
        mapped = _map_one(variant.position, position_map)
        if mapped is None:
            _warn_dropped(name, _describe(variant), "the position was removed")
            if variant.annot_id is not None:
                dropped_ids.add(variant.annot_id)
            continue
        simple.append(dataclasses.replace(variant, position=mapped))

    complexes: list[VariantComplex] = []
    for variant in entry.variant_complex:
        if pending[variant]:
            pending[variant] -= 1
            continue
        mapped_range = _remap_range(variant.start_pos, variant.end_pos, position_map)
        if mapped_range is None:
            _warn_dropped(name, _describe(variant), "a bound was removed")
            if variant.annot_id is not None:
                dropped_ids.add(variant.annot_id)
            continue
        complexes.append(dataclasses.replace(variant, start_pos=mapped_range[0], end_pos=mapped_range[1]))
    return tuple(simple), tuple(complexes)


def _resolvable(refs: Iterable[int | str], known_ids: set[int], gone_ids: set[int]) -> bool:
    """True unless a reference points at an annotation this form no longer has.

    A reference that was already dangling in the entry stays as it is: the reader warns
    about it, and this is not the place to repair the file.
    """
    for ref in refs:
        if isinstance(ref, int) and ref in gone_ids and ref in known_ids:
            return False
    return True


def _remap_refs(
    entry: SequenceEntry, known_ids: set[int], gone_ids: set[int], position_map: Sequence[int | None], name: str
) -> tuple[tuple[DisulfideBond, ...], tuple[Proteoform, ...]]:
    """Keep the annotations that reference ids (spec §3.4.2) only while they resolve."""
    bonds: list[DisulfideBond] = []
    for bond in entry.disulfide_bond:
        # annot_id_refs are annotation ids, not positions, so nothing is remapped.
        if not _resolvable(bond.annot_id_refs, known_ids, gone_ids):
            if bond.annot_id is not None:
                gone_ids.add(bond.annot_id)
            continue
        bonds.append(bond)

    proteoforms: list[Proteoform] = []
    for form in entry.proteoform:
        if not _resolvable(form.annot_id_refs, known_ids, gone_ids):
            continue
        ranges: list[SequenceRange] = []
        dropped = False
        for item in form.ranges:
            mapped = _remap_range(item.start, item.end, position_map)
            if mapped is None:
                _warn_dropped(
                    name, f"Proteoform {form.proteoform_id} range {item.start}-{item.end}", "a bound was removed"
                )
                dropped = True
                break
            ranges.append(SequenceRange(start=mapped[0], end=mapped[1]))
        if dropped:
            continue
        proteoforms.append(dataclasses.replace(form, ranges=tuple(ranges)))
    return tuple(bonds), tuple(proteoforms)


def _apply(entry: SequenceEntry, variants: Sequence[Variant]) -> VariantForm:
    """The body of :func:`entry_apply_variants`, without the ``errors`` handling."""
    name = f"{entry.prefix}:{entry.db_unique_id}"
    spans = _ordered(variants, len(entry.sequence), name)
    sequence, position_map, substituted, truncated_at = _splice(entry.sequence, spans)

    applied_counts = Counter(variant for _, _, variant in spans)
    known_ids = {
        annotation.annot_id
        for group in (
            entry.variant_simple,
            entry.variant_complex,
            entry.mod_res_unimod,
            entry.mod_res_psi,
            entry.mod_res,
            entry.processed,
            entry.disulfide_bond,
            entry.proteoform,
        )
        for annotation in group
        if annotation.annot_id is not None
    }
    # An applied variant is in the sequence, so it is no longer an annotation: a
    # \Proteoform that references it describes this form rather than one reachable from it.
    gone_ids = {variant.annot_id for _, _, variant in spans if variant.annot_id is not None}

    simple, complexes = _remap_pending_variants(entry, applied_counts, position_map, name, gone_ids)
    remapped = {
        "variant_simple": simple,
        "variant_complex": complexes,
        "mod_res_unimod": _remap_mods(entry.mod_res_unimod, position_map, substituted, name, gone_ids),
        "mod_res_psi": _remap_mods(entry.mod_res_psi, position_map, substituted, name, gone_ids),
        "mod_res": _remap_mods(entry.mod_res, position_map, substituted, name, gone_ids),
        "processed": _remap_processed(entry.processed, position_map, name, gone_ids),
    }
    bonds, proteoforms = _remap_refs(entry, known_ids, gone_ids, position_map, name)

    derived = dataclasses.replace(
        entry,
        sequence=sequence,
        length=None if entry.length is None else len(sequence),
        disulfide_bond=bonds,
        proteoform=proteoforms,
        custom_values=dict(entry.custom_values),
        extra=dict(entry.extra),
        **remapped,
    )
    return VariantForm(
        entry=derived,
        applied=tuple(variants),
        position_map=tuple(position_map),
        truncated_at=truncated_at,
        is_reference=not variants,
    )


def reference_form(entry: SequenceEntry) -> VariantForm:
    """The entry as written, with an identity position map."""
    return VariantForm(
        entry=entry,
        applied=(),
        position_map=(None, *range(1, len(entry.sequence) + 1)),
        truncated_at=None,
        is_reference=True,
    )


@overload
def entry_apply_variants(
    entry: SequenceEntry, variants: Iterable[Variant] = ..., *, errors: Literal["raise"] = ...
) -> VariantForm: ...


@overload
def entry_apply_variants(
    entry: SequenceEntry, variants: Iterable[Variant] = ..., *, errors: Literal["skip"]
) -> VariantForm | None: ...


def entry_apply_variants(
    entry: SequenceEntry,
    variants: Iterable[Variant] = (),
    *,
    errors: Literal["raise", "skip"] = "raise",
) -> VariantForm | None:
    """Apply *variants* to *entry*'s sequence. See ``SequenceEntry.apply_variants``."""
    if errors not in ("raise", "skip"):
        raise _error(f"unknown errors mode {errors!r}", "use errors='raise' or errors='skip'")
    chosen = tuple(variants)
    if not chosen:
        return reference_form(entry)
    try:
        return _apply(entry, chosen)
    except PeffError:
        if errors == "skip":
            return None
        raise


def _combination_count(total: int, max_combination: int, include_reference: bool) -> int:
    return (1 if include_reference else 0) + sum(math.comb(total, size) for size in range(1, max_combination + 1))


def entry_variant_forms(
    entry: SequenceEntry,
    *,
    max_combination: int = 1,
    proteoforms: bool = False,
    mandatory: bool = False,
    include_reference: bool = True,
    limit: int = DEFAULT_FORM_LIMIT,
) -> Iterator[VariantForm]:
    """Enumerate the variant forms of *entry*. See ``SequenceEntry.variant_forms``."""
    name = f"{entry.prefix}:{entry.db_unique_id}"
    variants: tuple[Variant, ...] = (*entry.variant_simple, *entry.variant_complex)

    if mandatory:
        # Spec section 3.4.1: in a ProteoformDb=true database every annotation MUST be
        # applied, so the entry describes exactly one sequence.
        yield _apply(entry, variants) if variants else reference_form(entry)
        return

    if proteoforms:
        # Spec section 3.4.2: the file lists the combinations it means, by annotation id.
        by_id = {variant.annot_id: variant for variant in variants if variant.annot_id is not None}
        for form in entry.proteoform:
            chosen = tuple(by_id[ref] for ref in form.annot_id_refs if ref in by_id)
            try:
                yield _apply(entry, chosen) if chosen else reference_form(entry)
            except PeffError as err:
                warnings.warn(f"proteoform {form.proteoform_id} skipped: {err}", PeffWarning, stacklevel=2)
        return

    if max_combination < 0:
        raise _error("max_combination must not be negative", "use 1 for one variant per form")
    size_cap = min(max_combination, len(variants))
    count = _combination_count(len(variants), size_cap, include_reference)
    if count > limit:
        # Counted before a single sequence is built, so this is immediate.
        raise _error(
            f"{name}: {len(variants)} variants give {count} forms at max_combination={max_combination}, "
            f"over the limit of {limit}",
            "raise limit=, lower max_combination=, or use proteoforms=True to enumerate only the "
            "combinations the file declares (spec section 3.4.2)",
        )

    if include_reference:
        yield reference_form(entry)
    for size in range(1, size_cap + 1):
        for chosen in itertools.combinations(variants, size):
            try:
                yield _apply(entry, chosen)
            except PeffError as err:
                warnings.warn(f"combination skipped: {err}", PeffWarning, stacklevel=2)
