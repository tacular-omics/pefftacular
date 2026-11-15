"""Frozen dataclass models for PEFF file structures."""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from datetime import date, time
from typing import Literal, overload

CustomFieldValue = str | int | float | bool | date | time

# ---------------------------------------------------------------------------
# Header types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CustomKeyDef:
    """Definition of a custom key in a PEFF database header."""

    key_name: str
    description: str
    concept_curie: str | None = None
    regexp: str | None = None
    field_names: tuple[str, ...] = ()
    field_types: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class OptionalTagDef:
    """Definition of an optional tag used in entry annotations."""

    tag: str
    description: str


@dataclass(frozen=True, slots=True)
class SpecificKeyDef:
    """A database-specific key declared in the header (PEFF 1.0 section 3.3.1).

    ``# SpecificKey=name:description`` declares the key and
    ``# SpecificValue=name:values`` gives its allowed values (for example
    ``(available|unsure|not available)``). Entries then use it as ``\\name=value``;
    those values are kept in ``SequenceEntry.extra``. ``description`` or ``value`` is
    None when the header has no ``SpecificKey`` or ``SpecificValue`` line for the key.
    """

    key_name: str
    description: str | None = None
    value: str | None = None


@dataclass(frozen=True, slots=True)
class DatabaseHeader:
    """Metadata block for a single database within a PEFF file."""

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    prefix: str | None = None
    db_name: str | None = None
    db_description: str | None = None
    db_version: str | None = None
    db_date: str | None = None
    db_sources: tuple[str, ...] = ()
    general_comments: tuple[str, ...] = ()
    number_of_entries: int | None = None
    sequence_type: str | None = None
    decoy: bool | None = None
    conversion: str | None = None
    has_annotation_identifiers: bool = False
    proteoform_db: bool = False
    custom_key_defs: tuple[CustomKeyDef, ...] = ()
    optional_tag_defs: tuple[OptionalTagDef, ...] = ()
    specific_keys: tuple[SpecificKeyDef, ...] = ()
    extra: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class FileHeader:
    """Top-level header for a PEFF file."""

    # Explicitly unhashable: its DatabaseHeaders hold an ``extra`` dict.
    __hash__ = None  # type: ignore[assignment]

    peff_version: str
    general_comments: tuple[str, ...] = ()
    databases: tuple[DatabaseHeader, ...] = ()


# ---------------------------------------------------------------------------
# Entry annotation types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VariantSimple:
    """Single amino-acid substitution."""

    position: int | str
    new_amino_acid: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class VariantComplex:
    """Multi-residue variant (insertion, deletion, or substitution)."""

    start_pos: int | str
    end_pos: int | str
    new_sequence: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModResUnimod:
    """Modification annotated with a UNIMOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModResPsi:
    """Modification annotated with a PSI-MOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModRes:
    """Modification with a generic accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class Processed:
    """Processed molecule annotation (signal peptide, mature protein, etc.)."""

    start_pos: int | str
    end_pos: int | str
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class DisulfideBond:
    """Disulfide bond linking two previously-declared modification annotations.

    Per spec section 3.4.2, a ``\\DisulfideBond`` references two prior
    ``ModResPsi`` (half-cystine) entries *by their annotation IDs* — not by
    sequence position. ``annot_id_refs`` therefore holds those referenced IDs.
    """

    annot_id_refs: tuple[int | str, ...]
    description: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class SequenceRange:
    """A contiguous start–end range within a sequence (1-based, inclusive)."""

    start: int | str
    end: int | str


@dataclass(frozen=True, slots=True)
class Proteoform:
    """A specific proteoform — a combination of processing events and modifications."""

    proteoform_id: str
    ranges: tuple[SequenceRange, ...]
    annot_id_refs: tuple[int, ...]
    name: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class CustomKeyValue:
    """A parsed value for a header-declared custom key on an entry."""

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    key_name: str
    fields: dict[str, CustomFieldValue] = field(default_factory=dict)
    raw: str = ""


# ---------------------------------------------------------------------------
# Sequence entry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SequenceEntry:
    """A single sequence entry in a PEFF file.

    Frozen and compared by value (``==``), but **not hashable**: ``custom_values`` and
    ``extra`` are dicts, so ``hash(entry)`` raises :class:`TypeError`. Key sets or dicts
    by ``(entry.prefix, entry.db_unique_id)`` instead.

    ``id`` holds the ``\\ID=`` key (spec §3.3.4); the name shadows the builtin but is kept
    for API stability.
    """

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    prefix: str
    db_unique_id: str
    sequence: str
    id: str | None = None
    db_unique_id_key: str | None = None
    pname: str | None = None
    gname: str | None = None
    ncbi_tax_id: int | None = None
    tax_name: str | None = None
    length: int | None = None
    sv: int | None = None
    ev: int | None = None
    pe: int | None = None
    decoy: bool | None = None
    comment: str | None = None
    variant_simple: tuple[VariantSimple, ...] = ()
    variant_complex: tuple[VariantComplex, ...] = ()
    mod_res_unimod: tuple[ModResUnimod, ...] = ()
    mod_res_psi: tuple[ModResPsi, ...] = ()
    mod_res: tuple[ModRes, ...] = ()
    processed: tuple[Processed, ...] = ()
    disulfide_bond: tuple[DisulfideBond, ...] = ()
    proteoform: tuple[Proteoform, ...] = ()
    custom_values: dict[str, tuple[CustomKeyValue, ...]] = field(default_factory=dict)
    extra: dict[str, str] = field(default_factory=dict)

    # -- conversions (logic in _convert.py; models stay plain data) ---------------------

    @classmethod
    def from_fasta(cls, header: str, sequence: str, *, prefix: str | None = None) -> "SequenceEntry":
        """Build an entry from a plain FASTA header and sequence.

        pefftacular has no dependencies, so this takes strings rather than a
        ``fastatacular.SequenceEntry``; pass ``(e.raw_header, e.sequence)`` for one of those.
        UniProt-style headers are split into fields: ``sp|P12345|NAME_HUMAN Name OS=Homo
        sapiens OX=9606 GN=ABC PE=1 SV=2`` gives ``prefix="sp"``, ``db_unique_id="P12345"``,
        ``id="NAME_HUMAN"``, ``pname``, ``tax_name``, ``ncbi_tax_id``, ``gname``, ``pe``,
        ``sv``; other ``KEY=value`` pairs go to ``extra``. ``length`` is set from the sequence.
        A repeated key keeps its first value, as in fastatacular, and the later occurrences
        are dropped.

        Args:
            header: The header line, with or without the leading ``>``.
            sequence: Residues; whitespace is removed.
            prefix: PEFF database prefix. Default: the ``db`` of a ``db:ACC`` (PEFF style,
                checked first) or ``db|ACC...`` identifier. Required for other identifiers.

        Raises:
            PeffError: The header is empty, or no prefix is given or derivable.
        """
        from pefftacular._convert import entry_from_fasta

        return entry_from_fasta(cls, header, sequence, prefix=prefix)

    def to_fasta(self) -> tuple[str, str]:
        """Return ``(header, sequence)`` for a plain FASTA record (header without ``>``).

        The header is ``prefix|db_unique_id[|id] [pname] [OS=] [OX=] [GN=] [PE=] [SV=]``
        followed by ``extra`` as ``KEY=value``. When ``db_unique_id`` contains ``|`` the
        identifier is the PEFF form ``prefix:db_unique_id`` instead, and ``id`` is dropped.
        Annotations (variants, modifications, processing, proteoforms), ``comment``, ``ev``
        and ``decoy`` have no FASTA form and are dropped.
        ``SequenceEntry.from_fasta(*entry.to_fasta())`` restores the fields above as long as
        ``pname`` and the values contain no `` KEY=`` text.
        """
        from pefftacular._convert import entry_to_fasta

        return entry_to_fasta(self)

    @overload
    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = ...,
        variants: Iterable[VariantSimple] = ...,
        errors: Literal["raise"] = ...,
    ) -> str: ...

    @overload
    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = ...,
        variants: Iterable[VariantSimple] = ...,
        errors: Literal["skip"],
    ) -> str | None: ...

    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = "psimod",
        variants: Iterable[VariantSimple] = (),
        errors: Literal["raise", "skip"] = "raise",
    ) -> str | None:
        """Render the sequence with its modifications as a ProForma 2.0 string.

        ``mods="psimod"`` writes ``\\ModResPsi`` sites (``S[MOD:00046]``) and
        ``"unimod"`` writes ``\\ModResUnimod`` sites (``S[UNIMOD:21]``). ``\\ModRes``
        sites are included when their accession is from the same vocabulary; others are
        left out. An empty accession is written by name (``[M:name]`` / ``[U:name]``). Every
        listed site is modified at once. Unknown positions (``?``) become a ProForma
        unknown-position prefix (``[MOD:00046]^2?SEQ``); they are dropped when a ``*``
        variant truncates the sequence, since they may lie in the removed part. A site
        listed in both ``\\ModResPsi``/``\\ModResUnimod`` and ``\\ModRes`` is written once, even in different case,
        spacing or zero-padding (``MOD:00046``, `` mod:00046``, ``MOD:46``, ``MOD: 00046``). The CV
        prefix is written in upper case, and numeric accessions in canonical form: PSI-MOD
        zero-padded to 5 digits (``MOD:00046``), Unimod unpadded (``UNIMOD:21``).
        PEFF cannot tell a terminal modification from one on the terminal residue, so all
        are written on the residue.

        Real files contain entries whose sites lie past the end of the sequence (12 of
        the 20,431 entries of the neXtProt human PEFF). Converting a whole file, pass
        ``errors="skip"`` to get ``None`` for those entries instead of an exception::

            forms = [p for e in entries if (p := e.to_proforma(errors="skip")) is not None]

        Args:
            mods: Which vocabulary to render.
            variants: ``VariantSimple`` substitutions to apply, usually a subset of
                ``self.variant_simple``. A modification on a substituted residue is
                dropped (spec section 3.3.10: a modified variant needs its own entry).
                ``*`` (stop) truncates the sequence before that position.
            errors: ``"raise"`` (default) raises :class:`PeffError` for an entry that
                cannot be written; ``"skip"`` returns ``None`` for it instead.

        Returns:
            The ProForma string, or ``None`` with ``errors="skip"`` when the entry
            cannot be written.

        Raises:
            PeffError: Unknown ``mods`` or ``errors`` (always). With ``errors="raise"``
                also a position outside the sequence, a non-numeric position other than
                ``?``, two different substitutions at one position, a residue ProForma
                cannot hold (``*``, ``-`` or anything outside ``A``-``Z``) in the written
                part of the sequence, or a square bracket, ``|`` or ``#`` in a modification
                accession or name (ProForma has no escape for them); the message
                starts with ``prefix:db_unique_id``.
        """
        from pefftacular._convert import entry_to_proforma

        return entry_to_proforma(self, mods=mods, variants=variants, errors=errors)

    def to_record(self) -> dict[str, str | int | bool | None]:
        """Return this entry as a flat ``dict`` (the keys of :func:`pefftacular.to_records`).

        Custom values are written from ``CustomKeyValue.raw`` (no header definitions are
        available here); :meth:`PeffReader.to_records` uses the file's definitions.
        """
        from pefftacular._records import entry_to_record

        return entry_to_record(self)

    # -- variants (logic in _variants.py; models stay plain data) -----------------------

    @overload
    def apply_variants(
        self,
        variants: Iterable["VariantSimple | VariantComplex"] = ...,
        *,
        errors: Literal["raise"] = ...,
    ) -> "VariantForm": ...

    @overload
    def apply_variants(
        self,
        variants: Iterable["VariantSimple | VariantComplex"] = ...,
        *,
        errors: Literal["skip"],
    ) -> "VariantForm | None": ...

    def apply_variants(
        self,
        variants: Iterable["VariantSimple | VariantComplex"] = (),
        *,
        errors: Literal["raise", "skip"] = "raise",
    ) -> "VariantForm | None":
        r"""Apply *variants* to the sequence and return the resulting :class:`VariantForm`.

        A PEFF entry's sequence line holds the reference sequence; ``\VariantSimple``
        (spec §3.3.8) and ``\VariantComplex`` (§3.3.9) describe sequences it does not
        contain. This builds one of them: the form's ``entry`` is an ordinary
        :class:`SequenceEntry`, so ``to_proforma``, ``to_record`` and ``write_peff`` work
        on it unchanged.

        Every other annotation is carried onto the new coordinates through
        ``VariantForm.position_map``:

        - A modification site the variant removed is dropped, as is one on a residue a
          ``VariantSimple`` substituted (§3.3.10: a modified variant needs its own entry).
          A ``?`` position is kept as written. A modification left with no site is dropped.
        - ``Processed``, un-applied ``VariantComplex`` and ``Proteoform`` ranges are
          remapped by both endpoints, and dropped rather than clipped when either endpoint
          was removed.
        - ``DisulfideBond`` and ``Proteoform`` references are annotation ids, not positions
          (§3.4.2), so they are copied; the annotation holding them is dropped when a
          reference no longer resolves, including a reference to a variant this form
          applied. A reference that was already dangling is left alone.
        - ``length`` is recomputed when the entry had one; the applied variants are
          removed from ``variant_simple`` / ``variant_complex``.

        Each drop emits a :class:`PeffWarning`; nothing is silent.

        Args:
            variants: The variants to apply, usually a subset of ``self.variant_simple``
                and ``self.variant_complex``. Empty gives the reference form.
            errors: ``"raise"`` (default) raises :class:`PeffError`; ``"skip"`` returns
                ``None`` instead.

        Returns:
            The :class:`VariantForm`, or ``None`` with ``errors="skip"``.

        Raises:
            PeffError: Unknown ``errors`` (always). With ``errors="raise"`` also a position
                that is not a number, a position outside the sequence, a ``VariantSimple``
                whose new residue is not a single letter or ``*``, or two variants that
                change the same position: the spec does not define how overlapping
                variants combine, so one per form.

        Example:
            >>> from pefftacular import SequenceEntry, VariantComplex, VariantSimple
            >>> entry = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MSTKW",
            ...                       variant_simple=(VariantSimple(2, "C"),),
            ...                       variant_complex=(VariantComplex(4, 5, ""),))
            >>> entry.apply_variants(entry.variant_simple).sequence
            'MCTKW'
            >>> entry.apply_variants(entry.variant_complex).position_map
            (None, 1, 2, 3, None, None)
            >>> tuple(form.sequence for form in entry.variant_forms())
            ('MSTKW', 'MCTKW', 'MST')
        """
        from pefftacular._variants import entry_apply_variants

        return entry_apply_variants(self, variants, errors=errors)

    def variant_forms(
        self,
        *,
        max_combination: int = 1,
        proteoforms: bool = False,
        mandatory: bool = False,
        include_reference: bool = True,
        limit: int = 10_000,
    ) -> "Iterator[VariantForm]":
        r"""Iterate over the sequence forms this entry describes.

        A variant is optional by default ("MAY be substituted", spec §3.3.8), so the
        default yields the reference form and then one form per variant - the linear
        reading, and the only one that stays bounded on real files. The spec's insulin
        example (§3.4.2, shipped as ``tests/fixtures/PEFF_AnnotID_Insulin_Valid.peff``)
        carries 70 ``VariantSimple`` annotations: 71 forms by default, 2,486 for every
        pair, and 2**70 for the full power set, against the 11 proteoforms the file
        itself declares.

        Args:
            max_combination: How many variants one form may apply. Combinations whose
                variants overlap are skipped with a :class:`PeffWarning`.
            proteoforms: Yield one form per ``\Proteoform`` definition instead, applying
                the variants it references by annotation id (§3.4.2). The reference form
                is not added.
            mandatory: Yield exactly one form with every variant applied, for an entry
                from a database whose header has ``ProteoformDb=true``, where annotations
                MUST all be applied (§3.4.1). Read the flag from
                ``DatabaseHeader.proteoform_db``.
            include_reference: Whether to yield the entry as written first.
            limit: Ceiling on the number of forms. The count is computed before any
                sequence is built, so an impossible request fails at once.

        Yields:
            One :class:`VariantForm` at a time; the iterator is lazy, so ``limit`` bounds
            work rather than memory.

        Raises:
            PeffError: ``max_combination`` is negative, or the requested number of forms
                is over ``limit``.
        """
        from pefftacular._variants import entry_variant_forms

        return entry_variant_forms(
            self,
            max_combination=max_combination,
            proteoforms=proteoforms,
            mandatory=mandatory,
            include_reference=include_reference,
            limit=limit,
        )


@dataclass(frozen=True, slots=True)
class VariantForm:
    """One sequence form of an entry, with its variants applied.

    Returned by :meth:`SequenceEntry.apply_variants` and
    :meth:`SequenceEntry.variant_forms`.

    Attributes:
        entry: The derived entry: the variant ``sequence``, the applied variants removed,
            every other annotation remapped. An ordinary :class:`SequenceEntry`.
        applied: The variants that were applied, in the order given.
        position_map: 1-based original position to 1-based position in
            ``entry.sequence``, ``None`` where the variant removed the residue. Index 0 is
            unused, so ``position_map[p]`` reads naturally.
        truncated_at: The 1-based position of a ``*`` (nonsense) substitution, before
            which the sequence ends (spec §3.3.8), or ``None``.
        is_reference: True for the entry as written, with no variant applied.
    """

    __hash__ = None  # type: ignore[assignment]

    entry: SequenceEntry
    applied: tuple[VariantSimple | VariantComplex, ...]
    position_map: tuple[int | None, ...]
    truncated_at: int | None = None
    is_reference: bool = False

    @property
    def sequence(self) -> str:
        """Shorthand for ``self.entry.sequence``."""
        return self.entry.sequence
