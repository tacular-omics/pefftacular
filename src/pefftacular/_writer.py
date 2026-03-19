"""PEFF file writer — serializes models back to PEFF format."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import IO

from pefftacular._models import (
    CustomKeyDef,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    Processed,
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


def _serialize_variant_simple(items: tuple[VariantSimple, ...]) -> str:
    parts: list[str] = []
    for v in items:
        fields = [_fmt_position(v.position), v.new_amino_acid]
        if v.tag:
            fields.append(v.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_variant_complex(items: tuple[VariantComplex, ...]) -> str:
    parts: list[str] = []
    for v in items:
        fields = [_fmt_position(v.start_pos), _fmt_position(v.end_pos), v.new_sequence]
        if v.tag:
            fields.append(v.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_mod_res_like(items: tuple[ModResUnimod | ModResPsi | ModRes, ...]) -> str:
    parts: list[str] = []
    for m in items:
        fields = [_fmt_positions(m.positions), m.accession, m.name]
        if m.tag:
            fields.append(m.tag)
        parts.append(f"({'|'.join(fields)})")
    return "".join(parts)


def _serialize_processed(items: tuple[Processed, ...]) -> str:
    parts: list[str] = []
    for p in items:
        fields = [_fmt_position(p.start_pos), _fmt_position(p.end_pos), p.accession, p.name]
        if p.tag:
            fields.append(p.tag)
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
        if db.db_version is not None:
            out.write(f"# DbVersion={db.db_version}\n")
        if db.db_source is not None:
            out.write(f"# DbSource={db.db_source}\n")
        if db.number_of_entries is not None:
            out.write(f"# NumberOfEntries={db.number_of_entries}\n")
        if db.sequence_type is not None:
            out.write(f"# SequenceType={db.sequence_type}\n")
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
            raise PeffWriteError(
                f"SequenceEntry has an empty prefix: db_unique_id={entry.db_unique_id!r}"
            )
        if not entry.db_unique_id:
            raise PeffWriteError(
                f"SequenceEntry has an empty db_unique_id: prefix={entry.prefix!r}"
            )
        if not entry.sequence:
            raise PeffWriteError(
                f"SequenceEntry {entry.prefix}:{entry.db_unique_id!r} has an empty sequence"
            )

    if isinstance(dest, (str, Path)):
        with Path(dest).open("w") as f:
            _write_header(header, f)
            for entry in entry_list:
                _write_entry(entry, f)
    else:
        _write_header(header, dest)
        for entry in entry_list:
            _write_entry(entry, dest)
