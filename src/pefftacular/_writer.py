"""PEFF file writer — serializes models back to PEFF format."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import IO

from pefftacular._models import (
    CustomKeyDef,
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
from pefftacular.errors import PeffWriteError

_SEQ_LINE_WIDTH = 60

# Canonical key order for entry description lines.
_ENTRY_KEY_ORDER = (
    "Length",
    "PName",
    "GName",
    "NcbiTaxId",
    "TaxName",
    "SV",
    "EV",
    "PE",
    "Decoy",
    "VariantSimple",
    "VariantComplex",
    "ModResUnimod",
    "ModResPsi",
    "ModRes",
    "Processed",
)


# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------


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
        fields = [_fmt_annot_position(v.annot_id, v.position), v.new_amino_acid]
        if v.tag:
            fields.append(v.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_variant_complex(items: tuple[VariantComplex, ...]) -> str:
    parts: list[str] = []
    for v in items:
        fields = [_fmt_annot_position(v.annot_id, v.start_pos), _fmt_position(v.end_pos), v.new_sequence]
        if v.tag:
            fields.append(v.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_mod_res_like(items: tuple[ModResUnimod | ModResPsi | ModRes, ...]) -> str:
    parts: list[str] = []
    for m in items:
        fields = [_fmt_annot_positions(m.annot_id, m.positions), m.accession, m.name]
        if m.tag:
            fields.append(m.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_processed(items: tuple[Processed, ...]) -> str:
    parts: list[str] = []
    for p in items:
        fields = [_fmt_annot_position(p.annot_id, p.start_pos), _fmt_position(p.end_pos), p.accession, p.name]
        if p.tag:
            fields.append(p.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_disulfide_bond(items: tuple[DisulfideBond, ...]) -> str:
    parts: list[str] = []
    for d in items:
        fields = [_fmt_annot_positions(d.annot_id, d.positions)]
        if d.description is not None:
            fields.append(d.description)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_proteoform(items: tuple[Proteoform, ...]) -> str:
    parts: list[str] = []
    for p in items:
        pf_id = f"{p.annot_id}:{p.proteoform_id}" if p.annot_id is not None else p.proteoform_id
        ranges_str = ",".join(f"{r.start}-{r.end}" for r in p.ranges)
        refs_str = ",".join(str(i) for i in p.annot_id_refs)
        fields = [pf_id, ranges_str, refs_str, p.name or ""]
        # strip trailing empty fields but keep at least 3
        while len(fields) > 3 and not fields[-1]:
            fields.pop()
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_custom_key_def(ckd: CustomKeyDef) -> str:
    """Serialize a single CustomKeyDef to its parenthesized form."""
    parts = [f"KeyName={ckd.key_name}", f'Description="{ckd.description}"']
    if ckd.regexp is not None:
        parts.append(f'RegExp="{ckd.regexp}"')
    if ckd.field_names:
        parts.append(f"FieldNames={','.join(ckd.field_names)}")
    if ckd.field_types:
        parts.append(f"FieldTypes={','.join(ckd.field_types)}")
    return f"({'|'.join(parts)})"


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


def _write_entry(entry: SequenceEntry, out: IO[str]) -> None:
    # Build key-value pairs in canonical order
    kv_parts: list[str] = []

    if entry.id is not None:
        kv_parts.append(f"\\ID={entry.id}")
    if entry.db_unique_id_key is not None:
        kv_parts.append(f"\\DbUniqueId={entry.db_unique_id_key}")
    if entry.length is not None:
        kv_parts.append(f"\\Length={entry.length}")
    if entry.pname is not None:
        kv_parts.append(f"\\PName={entry.pname}")
    if entry.gname is not None:
        kv_parts.append(f"\\GName={entry.gname}")
    if entry.ncbi_tax_id is not None:
        kv_parts.append(f"\\NcbiTaxId={entry.ncbi_tax_id}")
    if entry.tax_name is not None:
        kv_parts.append(f"\\TaxName={entry.tax_name}")
    if entry.sv is not None:
        kv_parts.append(f"\\SV={entry.sv}")
    if entry.ev is not None:
        kv_parts.append(f"\\EV={entry.ev}")
    if entry.pe is not None:
        kv_parts.append(f"\\PE={entry.pe}")
    if entry.decoy is not None:
        kv_parts.append(f"\\Decoy={'true' if entry.decoy else 'false'}")
    if entry.comment is not None:
        kv_parts.append(f"\\Comment={entry.comment}")
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
    for k, v in entry.extra.items():
        kv_parts.append(f"\\{k}={v}")

    desc_suffix = " ".join(kv_parts)
    if desc_suffix:
        out.write(f">{entry.prefix}:{entry.db_unique_id} {desc_suffix}\n")
    else:
        out.write(f">{entry.prefix}:{entry.db_unique_id}\n")

    # Sequence wrapped at 60 chars
    seq = entry.sequence
    for i in range(0, len(seq), _SEQ_LINE_WIDTH):
        out.write(seq[i : i + _SEQ_LINE_WIDTH] + "\n")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def write_peff(header: FileHeader, entries: Iterable[SequenceEntry], dest: str | Path | IO[str]) -> None:
    """Write a complete PEFF file."""
    if header is None:
        raise PeffWriteError("header must not be None")

    entry_list = list(entries)
    for entry in entry_list:
        if not entry.prefix:
            raise PeffWriteError(f"SequenceEntry has an empty prefix: db_unique_id={entry.db_unique_id!r}")
        if not entry.db_unique_id:
            raise PeffWriteError(f"SequenceEntry has an empty db_unique_id: prefix={entry.prefix!r}")
        if not entry.sequence:
            raise PeffWriteError(f"SequenceEntry {entry.prefix}:{entry.db_unique_id!r} has an empty sequence")

    if isinstance(dest, (str, Path)):
        with Path(dest).open("w") as f:
            _write_header(header, f)
            for entry in entry_list:
                _write_entry(entry, f)
    else:
        _write_header(header, dest)
        for entry in entry_list:
            _write_entry(entry, dest)
