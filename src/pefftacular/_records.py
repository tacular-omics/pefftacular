"""Flat ``dict`` records for building data frames (pandas, polars) without a dependency."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from pefftacular._compression import Compression
    from pefftacular._models import CustomKeyDef, FileHeader, SequenceEntry
    from pefftacular.errors import PeffWriteError

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
_RECORD_KEY_SET = frozenset(RECORD_KEYS)


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


def to_records(
    source: str | Path | IO[str] | IO[bytes] | Iterable[SequenceEntry], *, compression: Compression = "infer"
) -> list[Record]:
    """Return one flat ``dict`` per PEFF entry, ready for ``pandas.DataFrame(records)``.

    ``source`` is a path (``str`` or ``Path``, may be compressed), an open text handle,
    or an iterable of :class:`SequenceEntry` (e.g. the list from :func:`read_peff`).
    ``compression`` applies to a path or handle and is described on :class:`PeffReader`.
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
    from pefftacular._compression import check_compression
    from pefftacular._parser import PeffReader

    check_compression(compression)
    if isinstance(source, (str, Path)) or hasattr(source, "read"):
        with PeffReader(cast("str | Path | IO[str] | IO[bytes]", source), compression=compression) as reader:
            return reader.to_records()
    return [entry_to_record(e) for e in cast("Iterable[SequenceEntry]", source)]


# ---------------------------------------------------------------------------
# Records back to entries
# ---------------------------------------------------------------------------

# Column to PEFF key, in the order ``_writer._format_entry`` emits them, so the
# description line rebuilt below is the line ``write_peff`` would have written.
_SINGLE_KEYS: tuple[tuple[str, str], ...] = (
    ("id", "ID"),
    ("db_unique_id_key", "DbUniqueId"),
    ("length", "Length"),
    ("pname", "PName"),
    ("gname", "GName"),
    ("ncbi_tax_id", "NcbiTaxId"),
    ("tax_name", "TaxName"),
    ("sv", "SV"),
    ("ev", "EV"),
    ("pe", "PE"),
    ("decoy", "Decoy"),
    ("comment", "Comment"),
    ("variant_simple", "VariantSimple"),
    ("variant_complex", "VariantComplex"),
    ("mod_res_unimod", "ModResUnimod"),
    ("mod_res_psi", "ModResPsi"),
    ("mod_res", "ModRes"),
    ("processed", "Processed"),
    ("disulfide_bond", "DisulfideBond"),
    ("proteoform", "Proteoform"),
)
_IDENTITY_KEYS = ("prefix", "db_unique_id", "sequence")
_VERBATIM_KEYS = ("custom_values", "extra")
_INT_KEYS = frozenset({"ncbi_tax_id", "length", "sv", "ev", "pe"})
_ESCAPED_KEYS = frozenset({"pname", "gname", "tax_name", "comment"})
# Values the writer emits as given: an empty one would not read back as itself.
_NON_EMPTY_KEYS = frozenset(
    {
        "variant_simple",
        "variant_complex",
        "mod_res_unimod",
        "mod_res_psi",
        "mod_res",
        "processed",
        "disulfide_bond",
        "proteoform",
        *_IDENTITY_KEYS,
        *_VERBATIM_KEYS,
    }
)


def _absent(value: object) -> bool:
    """Is *value* a missing value? ``None`` and the frame markers are; ``""`` is not.

    ``""`` is a value a PEFF file can hold (``\\PName=`` parses back as ``""``), so it is
    kept, unlike in fastatacular where an empty string counts as absent.
    """
    if value is None:
        return True
    try:
        return bool(value != value)  # float NaN, numpy.nan, pandas.NaT
    except TypeError:  # pandas.NA: NA != NA is NA, whose truth value raises
        return True


def _rec_error(message: str, index: int, hint: str) -> PeffWriteError:
    from pefftacular.errors import PeffWriteError

    # PeffWriteError prefixes "Entry N: " from ``index``, which is the record's position here.
    return PeffWriteError(message, index=index, hint=hint)


def _rec_int(value: object, index: int, key: str) -> int:
    """An int column, also accepting an integral float (a frame column with gaps) or digits."""
    if isinstance(value, bool):
        raise _rec_error(f"{key}={value!r} is a bool", index, f"{key} is an integer")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            pass
    raise _rec_error(f"{key}={value!r} is not an integer", index, f'{key} takes an int, 9606.0 or "9606"')


def _rec_bool(value: object, index: int) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.strip().lower() in ("true", "false"):
        return value.strip().lower() == "true"
    raise _rec_error(f"decoy={value!r} is not a boolean", index, 'decoy takes True, False, 0, 1, "true" or "false"')


def _rec_text(value: object, index: int, key: str) -> str:
    if not isinstance(value, str):
        raise _rec_error(f"{key}={value!r} is not a string", index, f"{key} takes text")
    if not value and key in _NON_EMPTY_KEYS:
        raise _rec_error(f"{key} is empty", index, f"give {key} a value or None, which leaves the key out")
    return value


def _clean_record(record: Mapping[str, object], index: int) -> dict[str, Any]:
    """Every ``RECORD_KEYS`` column of *record*, converted, with a missing value as ``None``.

    An unknown key or a missing one is an error: a dropped column is how a caller silently
    loses data (``custom_values`` above all), which is the one thing this function exists
    to prevent.
    """
    if not isinstance(record, Mapping):
        raise _rec_error(f"{record!r} is not a mapping", index, 'pass df.to_dict("records"), one mapping per entry')
    unknown = [key for key in record if key not in _RECORD_KEY_SET]
    if unknown:
        raise _rec_error(
            f"unknown key(s) {', '.join(map(repr, unknown))}",
            index,
            "RECORD_KEYS lists every key a record may have; drop the extra column",
        )
    missing = [key for key in RECORD_KEYS if key not in record]
    if missing:
        raise _rec_error(
            f"missing key(s) {', '.join(map(repr, missing))}",
            index,
            "to_records() emits all of RECORD_KEYS; add the column with None rather than dropping it",
        )
    values: dict[str, Any] = {}
    for key in RECORD_KEYS:
        value = record[key]
        if _absent(value):
            if key in _IDENTITY_KEYS:
                raise _rec_error(f"{key} is missing", index, f"every entry needs a {key}")
            values[key] = None
        elif key in _INT_KEYS:
            values[key] = _rec_int(value, index, key)
        elif key == "decoy":
            values[key] = _rec_bool(value, index)
        else:
            values[key] = _rec_text(value, index, key)
    return values


def _description_line(values: Mapping[str, Any], index: int, defs: Mapping[str, CustomKeyDef]) -> str:
    """Rebuild the ``>Prefix:Id \\Key=value`` line for one cleaned record."""
    from pefftacular import _writer as w
    from pefftacular._lexer import split_description_keys

    parts: list[str] = []
    seen: dict[str, str] = {}
    for column, key in _SINGLE_KEYS:
        value = values[column]
        if value is None:
            continue
        if column in _ESCAPED_KEYS:
            text = w._escape_component(value)
        elif column == "decoy":
            text = "true" if value else "false"
        else:
            text = str(value)
        seen[key] = column
        parts.append(f"\\{key}={text}")
    for column in _VERBATIM_KEYS:
        text = values[column]
        if text is None:
            continue
        keys = split_description_keys(text)
        if not keys:
            raise _rec_error(
                f"{column}={text!r} holds no \\Key=value pair",
                index,
                f"{column} is the text to_records() emits, e.g. \\MyKey=(a|b)",
            )
        for key in keys:
            if key in seen:
                raise _rec_error(
                    f"{column} repeats \\{key}=, which column {seen[key]!r} already sets",
                    index,
                    "a key may be set once, by its own column or inside custom_values / extra",
                )
            if column == "custom_values" and key not in defs:
                raise _rec_error(
                    f"custom_values uses \\{key}=, which no CustomKeyDef declares",
                    index,
                    "pass the file header, from_records(records, header), so its CustomKeyDefs are known; "
                    "an undeclared key is read back as an unknown key, not a custom value",
                )
            seen[key] = column
        parts.append(text)
    suffix = " ".join(parts)
    identifier = f">{values['prefix']}:{values['db_unique_id']}"
    return f"{identifier} {suffix}" if suffix else identifier


def _record_entry(
    record: Mapping[str, object], index: int, defs_by_prefix: Mapping[str, Mapping[str, CustomKeyDef]]
) -> SequenceEntry:
    from pefftacular._parser import _parse_entry
    from pefftacular.errors import PeffError

    values = _clean_record(record, index)
    defs = defs_by_prefix.get(values["prefix"], {})
    description = _description_line(values, index, defs)
    try:
        entry = _parse_entry(description, [values["sequence"]], custom_key_defs=dict(defs))
    except PeffError as err:
        raise _rec_error(
            f"would not read back: {err}",
            index,
            "a value holds text the reader treats as syntax; fix the column named in the message",
        ) from err
    # The record the entry gives back must be the record that came in. Any column this
    # function failed to honour is an error here instead of data quietly lost.
    written = entry_to_record(entry, defs)
    if written != values:
        key = next(k for k in RECORD_KEYS if written[k] != values[k])
        raise _rec_error(
            f"{key}={values[key]!r} reads back as {written[key]!r}",
            index,
            "the text columns take the form to_records() writes; use the value in this message",
        )
    return entry


def from_records(records: Iterable[Mapping[str, object]], header: FileHeader | None = None) -> list[SequenceEntry]:
    """Build :class:`SequenceEntry` objects from flat records, the inverse of :func:`to_records`.

    Writes a data frame back to PEFF::

        header, entries = read_peff("in.peff")
        df = pandas.DataFrame(to_records(entries))
        write_peff(header, from_records(df.to_dict("records"), header), "out.peff")

    Each record is a mapping with **every** key in :data:`RECORD_KEYS`. A key that is not
    one of them, or one of them left out, raises: a dropped column is how a caller loses
    data without noticing. ``None``, float ``NaN`` and the pandas/polars missing markers
    mean the entry has no value for that key, so the key is left out of the entry; ``""``
    is an empty value and is kept, because a PEFF file can hold one (``\\PName=`` reads
    back as ``""``). ``prefix``, ``db_unique_id`` and ``sequence`` are required.

    ``header`` is only read for its :class:`CustomKeyDef` declarations, matched to each
    record by ``prefix``. It is needed whenever a record has ``custom_values``, because
    only a declared key is read back as a custom value rather than an unknown one, and
    passing it is what keeps those values from being silently dropped. The entries are
    returned on their own, as :func:`to_records` takes them, so the caller keeps owning
    the header and passes it to :func:`write_peff` unchanged.

    ``ncbi_tax_id``, ``length``, ``sv``, ``ev`` and ``pe`` accept ints, integral floats
    (``9606.0``, from a frame column with gaps) and digit strings. ``decoy`` accepts
    ``True``, ``0``, ``1`` and ``"true"``/``"false"``. The annotation columns,
    ``custom_values`` and ``extra`` are the PEFF value text :func:`to_records` emits, and
    are used as given.

    Every entry is parsed from the description line that :func:`write_peff` would have
    written, then turned back into a record and compared with the one that came in, so a
    column that could not be honoured raises instead of changing the data.

    Example:
        >>> from pefftacular import SequenceEntry, from_records, to_records
        >>> entry = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="ACDE", pname="Alpha")
        >>> record = to_records([entry])[0]
        >>> record.get("pname"), record.get("variant_simple")
        ('Alpha', None)
        >>> from_records([record]) == [entry]
        True

    Raises:
        PeffWriteError: A missing or unknown key, a bad value, a key set twice, a
            ``custom_values`` key no ``CustomKeyDef`` declares, or a value that would not
            read back as itself. ``index`` is the record's position.
    """
    # Built the way PeffReader builds it (_parser.py:1153), so a record is read back with
    # the same defs the file itself would have been read with.
    defs_by_prefix: dict[str, dict[str, CustomKeyDef]] = {}
    for db in header.databases if header is not None else ():
        if db.prefix and db.custom_key_defs:
            defs_by_prefix[db.prefix] = {ckd.key_name: ckd for ckd in db.custom_key_defs}
    return [_record_entry(record, i, defs_by_prefix) for i, record in enumerate(records)]
