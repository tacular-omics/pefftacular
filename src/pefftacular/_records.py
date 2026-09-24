"""Flat ``dict`` records for building data frames (pandas, polars) without a dependency."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from pefftacular._models import CustomKeyDef, SequenceEntry

Record = dict[str, str | int | bool | None]

# Column order of every record; the keys never change between entries or files.
RECORD_KEYS: tuple[str, ...] = (
    "prefix",
    "db_unique_id",
    "id",
    "db_unique_id_key",
    "pname",
    "gname",
    "ncbi_tax_id",
    "tax_name",
    "length",
    "sv",
    "ev",
    "pe",
    "decoy",
    "comment",
    "variant_simple",
    "variant_complex",
    "mod_res_unimod",
    "mod_res_psi",
    "mod_res",
    "processed",
    "disulfide_bond",
    "proteoform",
    "custom_values",
    "extra",
    "sequence",
)


def _text(items: tuple[Any, ...], fmt: Callable[[Any], str]) -> str | None:
    return fmt(items) if items else None


def entry_to_record(entry: SequenceEntry, defs: Mapping[str, CustomKeyDef] | None = None) -> Record:
    from pefftacular import _writer as w

    custom = (
        " ".join(
            f"\\{key}={w._serialize_custom_values(items, (defs or {}).get(key))}"
            for key, items in entry.custom_values.items()
        )
        or None
    )
    extra = " ".join(f"\\{k}={v}" for k, v in entry.extra.items()) or None
    return {
        "prefix": entry.prefix,
        "db_unique_id": entry.db_unique_id,
        "id": entry.id,
        "db_unique_id_key": entry.db_unique_id_key,
        "pname": entry.pname,
        "gname": entry.gname,
        "ncbi_tax_id": entry.ncbi_tax_id,
        "tax_name": entry.tax_name,
        "length": entry.length,
        "sv": entry.sv,
        "ev": entry.ev,
        "pe": entry.pe,
        "decoy": entry.decoy,
        "comment": entry.comment,
        "variant_simple": _text(entry.variant_simple, w._serialize_variant_simple),
        "variant_complex": _text(entry.variant_complex, w._serialize_variant_complex),
        "mod_res_unimod": _text(entry.mod_res_unimod, w._serialize_mod_res_like),
        "mod_res_psi": _text(entry.mod_res_psi, w._serialize_mod_res_like),
        "mod_res": _text(entry.mod_res, w._serialize_mod_res_like),
        "processed": _text(entry.processed, w._serialize_processed),
        "disulfide_bond": _text(entry.disulfide_bond, w._serialize_disulfide_bond),
        "proteoform": _text(entry.proteoform, w._serialize_proteoform),
        "custom_values": custom,
        "extra": extra,
        "sequence": entry.sequence,
    }


def to_records(source: str | Path | IO[str] | Iterable[SequenceEntry]) -> list[Record]:
    """Return one flat ``dict`` per PEFF entry, ready for ``pandas.DataFrame(records)``.

    ``source`` is a path (``str`` or ``Path``, may be compressed), an open text handle,
    or an iterable of :class:`SequenceEntry` (e.g. the list from :func:`read_peff`).
    The file header is not included. Every record has the same keys, in this order
    (``RECORD_KEYS``); a key whose value is absent is ``None``:

    - ``prefix``, ``db_unique_id``: the identifier ``prefix:db_unique_id`` (str).
    - ``id`` (``\\ID``), ``db_unique_id_key`` (``\\DbUniqueId``), ``pname``, ``gname``,
      ``tax_name``, ``comment``: str.
    - ``ncbi_tax_id``, ``length``, ``sv``, ``ev``, ``pe``: int. ``decoy``: bool.
    - ``variant_simple``, ``variant_complex``, ``mod_res_unimod``, ``mod_res_psi``,
      ``mod_res``, ``processed``, ``disulfide_bond``, ``proteoform``: the annotation's
      PEFF value text as the writer emits it, e.g. ``"(12|L)(30|*)"``, or ``None`` when
      the entry has none (so ``df["variant_simple"].notna()`` selects entries with
      variants).
    - ``custom_values``, ``extra``: header-declared custom keys and unknown keys as
      ``\\Key=value`` text separated by spaces.
    - ``sequence``: str.

    The package does not use or require pandas or polars; the records are plain dicts.
    """
    from pefftacular._parser import PeffReader

    if isinstance(source, (str, Path)) or hasattr(source, "read"):
        with PeffReader(cast("str | Path | IO[str]", source)) as reader:
            return reader.to_records()
    return [entry_to_record(e) for e in cast("Iterable[SequenceEntry]", source)]
