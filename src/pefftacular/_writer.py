"""PEFF file writer — serializes models back to PEFF format."""

from __future__ import annotations

import dataclasses
import logging
import warnings
from collections.abc import Iterable
from datetime import date, time
from pathlib import Path
from typing import IO

from pefftacular._models import (
    CustomKeyDef,
    CustomKeyValue,
    DisulfideBond,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    Processed,
    Proteoform,
    SequenceEntry,
    VariantComplex,
    VariantSimple,
)
from pefftacular._parser import _parse_custom_value
from pefftacular.errors import PeffError, PeffWriteError

logger = logging.getLogger("pefftacular.writer")

_SEQ_LINE_WIDTH = 60


# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------


def _escape_component(s: str) -> str:
    r"""Backslash-escape a free-text component for the entry description line.

    ``\`` and ``|`` are always escaped so they don't read as an escape lead-in
    or a component separator. Parentheses are escaped only when the component's
    parens are *unbalanced* — balanced pairs (e.g. ``N-linked (GlcNAc...)``) are
    left intact per spec section 3.3.3, which keeps common names readable.
    """
    escaped = s.replace("\\", "\\\\").replace("|", "\\|")
    depth = 0
    balanced = True
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                balanced = False
                break
    if depth != 0:
        balanced = False
    if not balanced:
        escaped = escaped.replace("(", "\\(").replace(")", "\\)")
    return escaped


def _fmt_position(p: int | str) -> str:
    return str(p)


def _fmt_positions(positions: tuple[int | str, ...]) -> str:
    return ",".join(_fmt_position(p) for p in positions)


def _fmt_annot_position(annot_id: int | None, p: int | str) -> str:
    return f"{annot_id}:{p}" if annot_id is not None else _fmt_position(p)


def _fmt_annot_positions(annot_id: int | None, positions: tuple[int | str, ...]) -> str:
    raw = _fmt_positions(positions)
    return f"{annot_id}:{raw}" if annot_id is not None else raw


def _serialize_variant_simple(items: tuple[VariantSimple, ...]) -> str:
    parts: list[str] = []
    for v in items:
        fields = [_fmt_annot_position(v.annot_id, v.position), _escape_component(v.new_amino_acid)]
        if v.tag:
            fields.append(_escape_component(v.tag))
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_variant_complex(items: tuple[VariantComplex, ...]) -> str:
    parts: list[str] = []
    for v in items:
        fields = [
            _fmt_annot_position(v.annot_id, v.start_pos),
            _fmt_position(v.end_pos),
            _escape_component(v.new_sequence),
        ]
        if v.tag:
            fields.append(_escape_component(v.tag))
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_mod_res_like(items: tuple[ModResUnimod | ModResPsi | ModRes, ...]) -> str:
    parts: list[str] = []
    for m in items:
        fields = [
            _fmt_annot_positions(m.annot_id, m.positions),
            _escape_component(m.accession),
            _escape_component(m.name),
        ]
        if m.tag:
            fields.append(_escape_component(m.tag))
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_processed(items: tuple[Processed, ...]) -> str:
    parts: list[str] = []
    for p in items:
        fields = [
            _fmt_annot_position(p.annot_id, p.start_pos),
            _fmt_position(p.end_pos),
            _escape_component(p.accession),
            _escape_component(p.name),
        ]
        if p.tag:
            fields.append(_escape_component(p.tag))
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_disulfide_bond(items: tuple[DisulfideBond, ...]) -> str:
    parts: list[str] = []
    for d in items:
        fields = [_fmt_annot_positions(d.annot_id, d.annot_id_refs)]
        if d.description:
            fields.append(_escape_component(d.description))
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_proteoform(items: tuple[Proteoform, ...]) -> str:
    parts: list[str] = []
    for p in items:
        escaped_id = _escape_component(p.proteoform_id)
        pf_id = f"{p.annot_id}:{escaped_id}" if p.annot_id is not None else escaped_id
        ranges_str = ",".join(f"{r.start}-{r.end}" for r in p.ranges)
        refs_str = ",".join(str(i) for i in p.annot_id_refs)
        fields = [pf_id, ranges_str, refs_str, _escape_component(p.name) if p.name else ""]
        # strip trailing empty fields but keep at least 3
        while len(fields) > 3 and not fields[-1]:
            fields.pop()
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _quote_keydef_value(value: str) -> str:
    """Escape ``"`` and ``\\`` inside a quoted CustomKeyDef value."""
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _serialize_custom_key_def(ckd: CustomKeyDef) -> str:
    """Serialize a single CustomKeyDef to its parenthesized form."""
    parts = [f"KeyName={ckd.key_name}", f'Description="{_quote_keydef_value(ckd.description)}"']
    if ckd.concept_curie is not None:
        parts.append(f"ConceptCURIE={ckd.concept_curie}")
    if ckd.regexp is not None:
        parts.append(f'RegExp="{_quote_keydef_value(ckd.regexp)}"')
    if ckd.field_names:
        parts.append(f"FieldNames={','.join(ckd.field_names)}")
    if ckd.field_types:
        parts.append(f"FieldTypes={','.join(ckd.field_types)}")
    return f"({'|'.join(parts)})"


def _fmt_custom_field(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (date, time)):
        return value.isoformat()
    return str(value)


def _raw_matches_fields(v: CustomKeyValue, ckd: CustomKeyDef | None) -> bool:
    """True if ``v.raw`` still encodes ``v.fields``, so it can be written verbatim.

    A parsed value keeps its original text in ``raw``; after
    ``dataclasses.replace(v, fields=...)`` that text is stale and the fields must be
    rebuilt. Without a def the raw text cannot be re-parsed, so it is trusted.
    """
    if ckd is None:
        return True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            reparsed = _parse_custom_value(f"({v.raw})", ckd)
        except PeffError:
            return False
    return len(reparsed) == 1 and reparsed[0].fields == v.fields


def _serialize_custom_values(items: tuple[CustomKeyValue, ...], ckd: CustomKeyDef | None) -> str:
    """Serialize the per-key tuple of CustomKeyValue back into its description form."""
    parts: list[str] = []
    for v in items:
        if v.raw and _raw_matches_fields(v, ckd):
            parts.append(f"({v.raw})")
            continue
        if ckd is not None and ckd.field_names:
            ordered = [_fmt_custom_field(v.fields[name]) for name in ckd.field_names if name in v.fields]
        else:
            ordered = [_fmt_custom_field(val) for val in v.fields.values()]
        # Without a RegExp the reader splits on unescaped '|' and unescapes each
        # field, so escape to match. A RegExp-controlled key sees the raw item.
        if ckd is None or ckd.regexp is None:
            ordered = [_escape_component(f) for f in ordered]
        item = "|".join(ordered)
        if (
            ckd is not None
            and ckd.regexp is not None
            and not _raw_matches_fields(CustomKeyValue(v.key_name, v.fields, raw=item), ckd)
        ):
            raise PeffWriteError(
                f"Custom key {ckd.key_name!r}: fields {v.fields!r} cannot be written through RegExp {ckd.regexp!r}",
                hint="The joined value must match the RegExp and keep its parentheses balanced; change the fields",
            )
        parts.append(f"({item})")
    return "".join(parts)


def _line_break_at(obj: object, path: str) -> str | None:
    """Return the path of the first string under *obj* holding a line break, else ``None``.

    Every PEFF value is written on a single line, so a ``\\n`` or ``\\r`` anywhere
    would silently split a header or description line.
    """
    if isinstance(obj, str):
        return path if "\n" in obj or "\r" in obj else None
    if isinstance(obj, dict):
        for k, v in obj.items():
            found = _line_break_at(k, f"{path}[{k!r}]") or _line_break_at(v, f"{path}[{k!r}]")
            if found:
                return found
    elif isinstance(obj, (tuple, list)):
        for i, v in enumerate(obj):
            found = _line_break_at(v, f"{path}[{i}]")
            if found:
                return found
    elif dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        for f in dataclasses.fields(obj):
            found = _line_break_at(getattr(obj, f.name), f"{path}.{f.name}")
            if found:
                return found
    return None


# ---------------------------------------------------------------------------
# Header writing
# ---------------------------------------------------------------------------


def _write_header(header: FileHeader, out: IO[str]) -> None:
    out.write(f"# PEFF {header.peff_version}\n")

    for comment in header.general_comments:
        out.write(f"# GeneralComment={comment}\n")

    for db in header.databases:
        out.write("# //\n")
        if db.db_name is not None:
            out.write(f"# DbName={db.db_name}\n")
        if db.prefix is not None:
            out.write(f"# Prefix={db.prefix}\n")
        if db.db_description is not None:
            out.write(f"# DbDescription={db.db_description}\n")
        for comment in db.general_comments:
            out.write(f"# GeneralComment={comment}\n")
        if db.db_version is not None:
            out.write(f"# DbVersion={db.db_version}\n")
        if db.db_date is not None:
            out.write(f"# DbDate={db.db_date}\n")
        for src in db.db_sources:
            out.write(f"# DbSource={src}\n")
        if db.number_of_entries is not None:
            out.write(f"# NumberOfEntries={db.number_of_entries}\n")
        if db.sequence_type is not None:
            out.write(f"# SequenceType={db.sequence_type}\n")
        if db.decoy is not None:
            out.write(f"# Decoy={'true' if db.decoy else 'false'}\n")
        if db.conversion is not None:
            out.write(f"# Conversion={db.conversion}\n")
        if db.has_annotation_identifiers:
            out.write("# HasAnnotationIdentifiers=true\n")
        if db.proteoform_db:
            out.write("# ProteoformDb=true\n")
        for otd in db.optional_tag_defs:
            out.write(f"# OptionalTagDef={otd.tag}:{otd.description}\n")
        for ckd in db.custom_key_defs:
            out.write(f"# CustomKeyDef={_serialize_custom_key_def(ckd)}\n")
        for k, v in db.extra.items():
            out.write(f"# {k}={v}\n")

    out.write("# //\n")


# ---------------------------------------------------------------------------
# Entry writing
# ---------------------------------------------------------------------------


def _format_entry(
    entry: SequenceEntry,
    defs_by_prefix: dict[str, dict[str, CustomKeyDef]] | None = None,
) -> str:
    """Serialize one entry (description line plus wrapped sequence) to PEFF text."""
    # Build key-value pairs in canonical order
    kv_parts: list[str] = []

    if entry.id is not None:
        kv_parts.append(f"\\ID={entry.id}")
    if entry.db_unique_id_key is not None:
        kv_parts.append(f"\\DbUniqueId={entry.db_unique_id_key}")
    if entry.length is not None:
        kv_parts.append(f"\\Length={entry.length}")
    if entry.pname is not None:
        kv_parts.append(f"\\PName={_escape_component(entry.pname)}")
    if entry.gname is not None:
        kv_parts.append(f"\\GName={_escape_component(entry.gname)}")
    if entry.ncbi_tax_id is not None:
        kv_parts.append(f"\\NcbiTaxId={entry.ncbi_tax_id}")
    if entry.tax_name is not None:
        kv_parts.append(f"\\TaxName={_escape_component(entry.tax_name)}")
    if entry.sv is not None:
        kv_parts.append(f"\\SV={entry.sv}")
    if entry.ev is not None:
        kv_parts.append(f"\\EV={entry.ev}")
    if entry.pe is not None:
        kv_parts.append(f"\\PE={entry.pe}")
    if entry.decoy is not None:
        kv_parts.append(f"\\Decoy={'true' if entry.decoy else 'false'}")
    if entry.comment is not None:
        kv_parts.append(f"\\Comment={_escape_component(entry.comment)}")
    if entry.variant_simple:
        kv_parts.append(f"\\VariantSimple={_serialize_variant_simple(entry.variant_simple)}")
    if entry.variant_complex:
        kv_parts.append(f"\\VariantComplex={_serialize_variant_complex(entry.variant_complex)}")
    if entry.mod_res_unimod:
        kv_parts.append(f"\\ModResUnimod={_serialize_mod_res_like(entry.mod_res_unimod)}")
    if entry.mod_res_psi:
        kv_parts.append(f"\\ModResPsi={_serialize_mod_res_like(entry.mod_res_psi)}")
    if entry.mod_res:
        kv_parts.append(f"\\ModRes={_serialize_mod_res_like(entry.mod_res)}")
    if entry.processed:
        kv_parts.append(f"\\Processed={_serialize_processed(entry.processed)}")
    if entry.disulfide_bond:
        kv_parts.append(f"\\DisulfideBond={_serialize_disulfide_bond(entry.disulfide_bond)}")
    if entry.proteoform:
        kv_parts.append(f"\\Proteoform={_serialize_proteoform(entry.proteoform)}")
    if entry.custom_values:
        defs = (defs_by_prefix or {}).get(entry.prefix, {})
        for key, items in entry.custom_values.items():
            kv_parts.append(f"\\{key}={_serialize_custom_values(items, defs.get(key))}")
    for k, v in entry.extra.items():
        kv_parts.append(f"\\{k}={v}")

    desc_suffix = " ".join(kv_parts)
    lines = [
        f">{entry.prefix}:{entry.db_unique_id} {desc_suffix}"
        if desc_suffix
        else f">{entry.prefix}:{entry.db_unique_id}"
    ]

    # Sequence wrapped at 60 chars
    seq = entry.sequence
    lines.extend(seq[i : i + _SEQ_LINE_WIDTH] for i in range(0, len(seq), _SEQ_LINE_WIDTH))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def write_peff(header: FileHeader, entries: Iterable[SequenceEntry], dest: str | Path | IO[str]) -> None:
    """Write a complete PEFF file."""
    if header is None:
        raise PeffWriteError("header must not be None", hint="Pass a FileHeader instance, e.g. from read_peff()")

    bad = _line_break_at(header, "header")
    if bad:
        raise PeffWriteError(
            f"{bad} contains a line break", hint="PEFF header values are single-line; remove the \\n or \\r"
        )

    entry_list = list(entries)
    defs_by_prefix: dict[str, dict[str, CustomKeyDef]] = {}
    for db in header.databases:
        if db.prefix and db.custom_key_defs:
            defs_by_prefix[db.prefix] = {ckd.key_name: ckd for ckd in db.custom_key_defs}

    texts: list[str] = []
    for index, entry in enumerate(entry_list):
        if not entry.prefix:
            raise PeffWriteError(
                f"SequenceEntry has an empty prefix: db_unique_id={entry.db_unique_id!r}",
                index=index,
                hint="Every SequenceEntry needs a non-empty prefix matching a database in the header",
            )
        if not entry.db_unique_id:
            raise PeffWriteError(
                f"SequenceEntry has an empty db_unique_id: prefix={entry.prefix!r}",
                index=index,
                hint="Every SequenceEntry needs a non-empty db_unique_id (the accession after the prefix)",
            )
        if not entry.sequence:
            raise PeffWriteError(
                f"SequenceEntry {entry.prefix}:{entry.db_unique_id!r} has an empty sequence",
                index=index,
                hint="A PEFF entry must carry at least one residue in its sequence",
            )
        if ":" in entry.prefix or any(c.isspace() for c in entry.prefix):
            raise PeffWriteError(
                f"SequenceEntry prefix {entry.prefix!r} contains ':' or whitespace",
                index=index,
                hint="The prefix is the token before the first ':' of '>prefix:DbUniqueId'",
            )
        if any(c.isspace() for c in entry.db_unique_id):
            raise PeffWriteError(
                f"SequenceEntry db_unique_id {entry.db_unique_id!r} contains whitespace",
                index=index,
                hint="The DbUniqueId ends at the first space of the description line",
            )
        if ">" in entry.sequence or any(c.isspace() for c in entry.sequence):
            raise PeffWriteError(
                f"SequenceEntry {entry.prefix}:{entry.db_unique_id} sequence contains whitespace or '>'",
                index=index,
                hint="Pass the residues only; the writer wraps the sequence itself",
            )
        bad = _line_break_at(entry, "SequenceEntry")
        if bad:
            raise PeffWriteError(
                f"{entry.prefix}:{entry.db_unique_id}: {bad} contains a line break",
                index=index,
                hint="A PEFF description line is single-line; remove the \\n or \\r",
            )
        # Serialize now so a late failure (e.g. a RegExp-controlled custom key) is
        # raised before anything is written.
        try:
            texts.append(_format_entry(entry, defs_by_prefix))
        except PeffWriteError as err:
            raise PeffWriteError(str(err), index=index, hint=err.hint) from err

    if isinstance(dest, (str, Path)):
        logger.debug("writing PEFF file: %s (%d entries)", dest, len(entry_list))
        with Path(dest).open("w", encoding="utf-8") as f:
            _write_header(header, f)
            f.writelines(texts)
    else:
        logger.debug("writing PEFF to in-memory stream: %s (%d entries)", type(dest).__name__, len(entry_list))
        _write_header(header, dest)
        dest.writelines(texts)

    logger.info("write_peff: wrote %d entries across %d database(s)", len(entry_list), len(header.databases))
