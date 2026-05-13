"""PEFF file parser — header and sequence entry reading."""

from __future__ import annotations

import itertools
import re
import warnings
from collections.abc import Iterable, Iterator, Mapping
from datetime import date, time
from pathlib import Path
from typing import IO, Self

from pefftacular._lexer import split_description_keys, split_fields, split_items
from pefftacular._models import (
    CustomFieldValue,
    CustomKeyDef,
    CustomKeyValue,
    DatabaseHeader,
    DisulfideBond,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    OptionalTagDef,
    Processed,
    Proteoform,
    SequenceEntry,
    SequenceRange,
    VariantComplex,
    VariantSimple,
)
from pefftacular.errors import PeffParseError

# Compiled regex cache for CustomKeyDef.regexp.
_REGEXP_CACHE: dict[str, re.Pattern[str]] = {}
_ENUM_RE = re.compile(r"^enumeration\((.*)\)$")

# ---------------------------------------------------------------------------
# Position helpers
# ---------------------------------------------------------------------------

_HEADER_RE = re.compile(r"^# PEFF (\d+\.\d+)\s*$")


def parse_position(s: str) -> int | str:
    """Convert a position string to int, or leave as str if non-numeric (e.g. '?')."""
    try:
        return int(s)
    except ValueError:
        return s


def parse_positions(s: str) -> tuple[int | str, ...]:
    """Parse a comma-separated list of positions."""
    return tuple(parse_position(p.strip()) for p in s.split(","))


def _extract_annot_id(s: str) -> tuple[int | None, str]:
    """Split 'annotID:rest' into (annotID, rest). Returns (None, s) if no valid prefix."""
    colon = s.find(":")
    if colon > 0:
        try:
            return int(s[:colon]), s[colon + 1 :]
        except ValueError:
            pass
    return None, s


# ---------------------------------------------------------------------------
# Annotation parsers
# ---------------------------------------------------------------------------


def _parse_variant_simple(raw: str) -> tuple[VariantSimple, ...]:
    items = split_items(raw)
    results: list[VariantSimple] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 2:
            raise PeffParseError(f"VariantSimple needs >= 2 fields, got {len(fields)}", context=item)
        annot_id, pos_str = _extract_annot_id(fields[0])
        tag = fields[2] if len(fields) > 2 and fields[2] else None
        results.append(
            VariantSimple(position=parse_position(pos_str), new_amino_acid=fields[1], tag=tag, annot_id=annot_id)
        )
    return tuple(results)


def _parse_variant_complex(raw: str) -> tuple[VariantComplex, ...]:
    items = split_items(raw)
    results: list[VariantComplex] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 3:
            raise PeffParseError(f"VariantComplex needs >= 3 fields, got {len(fields)}", context=item)
        annot_id, pos_str = _extract_annot_id(fields[0])
        tag = fields[3] if len(fields) > 3 and fields[3] else None
        results.append(
            VariantComplex(
                start_pos=parse_position(pos_str),
                end_pos=parse_position(fields[1]),
                new_sequence=fields[2],
                tag=tag,
                annot_id=annot_id,
            )
        )
    return tuple(results)


def _parse_mod_res_like(raw: str) -> list[tuple[tuple[int | str, ...], str, str, str | None, int | None]]:
    """Shared parser for ModResUnimod / ModResPsi / ModRes."""
    items = split_items(raw)
    results: list[tuple[tuple[int | str, ...], str, str, str | None, int | None]] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 3:
            raise PeffParseError(f"ModRes-like key needs >= 3 fields, got {len(fields)}", context=item)
        annot_id, positions_str = _extract_annot_id(fields[0])
        positions = parse_positions(positions_str)
        accession = fields[1]
        name = fields[2]
        tag = fields[3] if len(fields) > 3 and fields[3] else None
        results.append((positions, accession, name, tag, annot_id))
    return results


def _parse_mod_res_unimod(raw: str) -> tuple[ModResUnimod, ...]:
    return tuple(
        ModResUnimod(positions=p, accession=a, name=n, tag=t, annot_id=i) for p, a, n, t, i in _parse_mod_res_like(raw)
    )


def _parse_mod_res_psi(raw: str) -> tuple[ModResPsi, ...]:
    return tuple(
        ModResPsi(positions=p, accession=a, name=n, tag=t, annot_id=i) for p, a, n, t, i in _parse_mod_res_like(raw)
    )


def _parse_mod_res(raw: str) -> tuple[ModRes, ...]:
    return tuple(
        ModRes(positions=p, accession=a, name=n, tag=t, annot_id=i) for p, a, n, t, i in _parse_mod_res_like(raw)
    )


def _parse_processed(raw: str) -> tuple[Processed, ...]:
    items = split_items(raw)
    results: list[Processed] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 4:
            raise PeffParseError(f"Processed needs >= 4 fields, got {len(fields)}", context=item)
        annot_id, pos_str = _extract_annot_id(fields[0])
        tag = fields[4] if len(fields) > 4 and fields[4] else None
        results.append(
            Processed(
                start_pos=parse_position(pos_str),
                end_pos=parse_position(fields[1]),
                accession=fields[2],
                name=fields[3],
                tag=tag,
                annot_id=annot_id,
            )
        )
    return tuple(results)


def _parse_disulfide_bond(raw: str) -> tuple[DisulfideBond, ...]:
    items = split_items(raw)
    results: list[DisulfideBond] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 1:
            raise PeffParseError("DisulfideBond needs >= 1 field", context=item)
        annot_id, positions_str = _extract_annot_id(fields[0])
        positions = parse_positions(positions_str)
        description = fields[1] if len(fields) > 1 and fields[1] else None
        results.append(DisulfideBond(positions=positions, description=description, annot_id=annot_id))
    return tuple(results)


def _parse_sequence_range(s: str) -> SequenceRange:
    """Parse 'start-end' into a SequenceRange."""
    if "-" in s:
        parts = s.split("-", 1)
        return SequenceRange(start=parse_position(parts[0]), end=parse_position(parts[1]))
    # fallback: treat whole string as start with unknown end
    return SequenceRange(start=parse_position(s), end=parse_position(s))


def _parse_proteoform(raw: str) -> tuple[Proteoform, ...]:
    items = split_items(raw)
    results: list[Proteoform] = []
    for item in items:
        fields = split_fields(item)
        if len(fields) < 2:
            raise PeffParseError("Proteoform needs >= 2 fields", context=item)
        annot_id, pf_id = _extract_annot_id(fields[0])
        # field[1]: comma-separated ranges like "1-110" or "90-110,25-54"
        range_strs = [r.strip() for r in fields[1].split(",") if r.strip()]
        ranges = tuple(_parse_sequence_range(r) for r in range_strs)
        # field[2]: comma-separated annotation ID refs (may be empty)
        annot_id_refs: tuple[int, ...] = ()
        if len(fields) > 2 and fields[2]:
            try:
                annot_id_refs = tuple(int(x.strip()) for x in fields[2].split(",") if x.strip())
            except ValueError as err:
                raise PeffParseError(f"Invalid annotation ID refs in Proteoform: {fields[2]!r}", context=item) from err
        name = fields[3] if len(fields) > 3 and fields[3] else None
        results.append(
            Proteoform(proteoform_id=pf_id, ranges=ranges, annot_id_refs=annot_id_refs, name=name, annot_id=annot_id)
        )
    return tuple(results)


# ---------------------------------------------------------------------------
# CustomKeyDef parser
# ---------------------------------------------------------------------------


def _parse_keydef_field(fr: str) -> tuple[str, str] | None:
    """Parse one ``defKey=value`` or ``defKey="value"`` token from a CustomKeyDef item.

    The value may be quoted with ``"..."`` (with ``\\"`` escaping inner quotes) or
    bare. Returns ``None`` for malformed fragments rather than raising.
    """
    s = fr.strip()
    eq = s.find("=")
    if eq <= 0:
        return None
    key = s[:eq]
    rest = s[eq + 1 :]
    if rest.startswith('"'):
        # Walk until matching close-quote. Only ``\"`` and ``\\`` are treated as
        # escapes; other backslashes are passed through so embedded regexes keep
        # their meaning.
        out: list[str] = []
        i = 1
        while i < len(rest):
            ch = rest[i]
            if ch == "\\" and i + 1 < len(rest) and rest[i + 1] in ('"', "\\"):
                out.append(rest[i + 1])
                i += 2
                continue
            if ch == '"':
                break
            out.append(ch)
            i += 1
        return key, "".join(out)
    return key, rest


def _parse_custom_key_def(raw: str) -> tuple[CustomKeyDef, ...]:
    """Parse a ``# CustomKeyDef=(KeyName=X|Description="Y"|...)`` line value.

    Multiple parenthesized items on the same line are supported, though the
    spec example uses one item per line.
    """
    items = split_items(raw)
    results: list[CustomKeyDef] = []
    for item in items:
        fields_raw = split_fields(item)
        kv: dict[str, str] = {}
        for fr in fields_raw:
            parsed = _parse_keydef_field(fr)
            if parsed is not None:
                kv[parsed[0]] = parsed[1]
        fn = tuple(kv["FieldNames"].split(",")) if kv.get("FieldNames") else ()
        ft = tuple(kv["FieldTypes"].split(",")) if kv.get("FieldTypes") else ()
        results.append(
            CustomKeyDef(
                key_name=kv.get("KeyName", ""),
                description=kv.get("Description", ""),
                concept_curie=kv.get("ConceptCURIE"),
                regexp=kv.get("RegExp"),
                field_names=fn,
                field_types=ft,
            )
        )
    return tuple(results)


# ---------------------------------------------------------------------------
# Custom-key value parser (entry side)
# ---------------------------------------------------------------------------


def _coerce_field(raw: str, xsd_type: str | None, *, key_name: str, field_name: str) -> CustomFieldValue:
    """Coerce a captured field string per the declared XSD type.

    On coercion failure, emit a warning and return the raw string. ``None``
    type or unknown type falls back to string.
    """
    if xsd_type is None:
        return raw
    t = xsd_type.strip()
    try:
        if t == "integer":
            return int(raw)
        if t == "decimal":
            return float(raw)
        if t == "boolean":
            return raw.strip().lower() in ("true", "1", "yes")
        if t == "date":
            return date.fromisoformat(raw)
        if t == "time":
            return time.fromisoformat(raw)
        if t == "string":
            return raw
        m = _ENUM_RE.match(t)
        if m:
            allowed = {opt.strip() for opt in m.group(1).split("|") if opt.strip()}
            if raw not in allowed:
                warnings.warn(
                    f"Custom key {key_name!r}: field {field_name!r}={raw!r} not in enumeration {sorted(allowed)}",
                    UserWarning,
                    stacklevel=3,
                )
            return raw
    except ValueError:
        warnings.warn(
            f"Custom key {key_name!r}: cannot coerce field {field_name!r}={raw!r} to {t}",
            UserWarning,
            stacklevel=3,
        )
        return raw
    return raw


def _compile_regexp(pattern: str) -> re.Pattern[str] | None:
    """Compile and cache a CustomKeyDef regexp. Returns None on bad pattern."""
    cached = _REGEXP_CACHE.get(pattern)
    if cached is not None:
        return cached
    try:
        compiled = re.compile(pattern)
    except re.error as err:
        warnings.warn(f"Invalid CustomKeyDef RegExp {pattern!r}: {err}", UserWarning, stacklevel=3)
        return None
    _REGEXP_CACHE[pattern] = compiled
    return compiled


def _parse_custom_value(raw: str, ckd: CustomKeyDef) -> tuple[CustomKeyValue, ...]:
    """Parse the value of a declared custom key into one or more CustomKeyValues.

    The value may be a single bare token or a sequence of parenthesized items.
    If the def has a ``RegExp``, it's applied per item and groups are mapped
    onto ``FieldNames``. Otherwise the item is split on ``|`` and zipped with
    ``FieldNames``. ``FieldTypes`` drives coercion in either case.
    """
    items = split_items(raw)
    results: list[CustomKeyValue] = []

    pattern = _compile_regexp(ckd.regexp) if ckd.regexp else None

    for item in items:
        fields: dict[str, CustomFieldValue] = {}
        if pattern is not None:
            m = pattern.fullmatch(item)
            if m is None:
                warnings.warn(
                    f"Custom key {ckd.key_name!r}: value {item!r} does not match RegExp {ckd.regexp!r}",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                groups = m.groups()
                for idx, raw_group in enumerate(groups):
                    if raw_group is None:
                        continue
                    fname = ckd.field_names[idx] if idx < len(ckd.field_names) else f"field{idx + 1}"
                    ftype = ckd.field_types[idx] if idx < len(ckd.field_types) else None
                    fields[fname] = _coerce_field(raw_group, ftype, key_name=ckd.key_name, field_name=fname)
        else:
            parts = split_fields(item)
            for idx, raw_part in enumerate(parts):
                if idx >= len(ckd.field_names) and not ckd.field_names:
                    fname = f"field{idx + 1}"
                elif idx < len(ckd.field_names):
                    fname = ckd.field_names[idx]
                else:
                    fname = f"field{idx + 1}"
                ftype = ckd.field_types[idx] if idx < len(ckd.field_types) else None
                fields[fname] = _coerce_field(raw_part, ftype, key_name=ckd.key_name, field_name=fname)

        results.append(CustomKeyValue(key_name=ckd.key_name, fields=fields, raw=item))

    return tuple(results)


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

_MANDATORY_DB_KEYS = {"Prefix", "DbVersion", "DbSource", "NumberOfEntries", "SequenceType"}


def _parse_file_header(lines: Iterable[str]) -> tuple[FileHeader, Iterator[str]]:
    """Parse the header section and return (FileHeader, remaining_lines_iterator)."""
    line_iter = iter(lines)
    line_no = 0

    # --- PEFF version line ---
    first_line = ""
    for raw_line in line_iter:
        line_no += 1
        first_line = raw_line.rstrip("\n\r")
        if first_line.strip():
            break
    else:
        raise PeffParseError("Empty file — no PEFF header found")

    m = _HEADER_RE.match(first_line)
    if not m:
        raise PeffParseError("First line must be '# PEFF <version>'", line=line_no, context=first_line)

    peff_version = m.group(1)
    if peff_version != "1.0":
        warnings.warn(f"PEFF version {peff_version} detected; only 1.0 is fully supported", stacklevel=2)

    # --- General comments and database blocks ---
    general_comments: list[str] = []
    databases: list[DatabaseHeader] = []
    current_db_lines: list[tuple[int, str]] = []
    in_general = True  # before the first "//"

    for raw_line in line_iter:
        line_no += 1
        line = raw_line.rstrip("\n\r")

        if not line.startswith("# ") and line != "#":
            # End of header — push this line back
            remaining = itertools.chain([raw_line], line_iter)
            break

        content = line[2:]  # strip "# "

        if content == "//":
            if in_general:
                in_general = False
            else:
                # Close the current database block
                databases.append(_build_database_header(current_db_lines))
                current_db_lines = []
            continue

        if in_general:
            if content.startswith("GeneralComment="):
                general_comments.append(content[len("GeneralComment=") :])
            # Other pre-// comment lines are ignored
        else:
            current_db_lines.append((line_no, content))
    else:
        # File ended during header
        if current_db_lines:
            databases.append(_build_database_header(current_db_lines))
        remaining = iter([])

    header = FileHeader(
        peff_version=peff_version,
        general_comments=tuple(general_comments),
        databases=tuple(databases),
    )
    return header, remaining


def _build_database_header(lines: list[tuple[int, str]]) -> DatabaseHeader:
    """Build a DatabaseHeader from the key=value lines inside a ``# //`` block."""
    # Collect values; multi-value keys accumulate into lists
    multi: dict[str, list[str]] = {}
    single: dict[str, str] = {}
    for _line_no, content in lines:
        eq_idx = content.find("=")
        if eq_idx == -1:
            continue
        key = content[:eq_idx]
        value = content[eq_idx + 1 :]
        if key in ("DbSource", "OptionalTagDef", "CustomKeyDef", "GeneralComment", "SpecificKey", "SpecificValue"):
            multi.setdefault(key, []).append(value)
        else:
            single[key] = value

    # Warn on missing mandatory keys
    all_keys = set(single) | set(multi)
    for mk in _MANDATORY_DB_KEYS:
        if mk not in all_keys:
            warnings.warn(f"Database header missing mandatory key: {mk}", stacklevel=3)

    custom_key_defs: tuple[CustomKeyDef, ...] = ()
    if "CustomKeyDef" in multi:
        accumulated: list[CustomKeyDef] = []
        for raw in multi.pop("CustomKeyDef"):
            accumulated.extend(_parse_custom_key_def(raw))
        custom_key_defs = tuple(accumulated)

    optional_tag_defs: tuple[OptionalTagDef, ...] = ()
    if "OptionalTagDef" in multi:
        defs: list[OptionalTagDef] = []
        for raw in multi.pop("OptionalTagDef"):
            colon = raw.find(":")
            if colon > 0:
                defs.append(OptionalTagDef(tag=raw[:colon], description=raw[colon + 1 :]))
            else:
                defs.append(OptionalTagDef(tag=raw, description=""))
        optional_tag_defs = tuple(defs)

    db_sources = tuple(multi.pop("DbSource", []))

    # Remove known multi keys that we don't expose further
    for k in ("GeneralComment", "SpecificKey", "SpecificValue"):
        multi.pop(k, None)

    prefix = single.pop("Prefix", None)
    db_name = single.pop("DbName", None)
    db_description = single.pop("DbDescription", None)
    db_version = single.pop("DbVersion", None)
    db_date = single.pop("DbDate", None)
    conversion = single.pop("Conversion", None)

    num_str = single.pop("NumberOfEntries", None)
    number_of_entries: int | None = None
    if num_str is not None:
        try:
            number_of_entries = int(num_str)
        except ValueError as err:
            raise PeffParseError(f"Invalid integer for NumberOfEntries: {num_str!r}", context=num_str) from err

    sequence_type = single.pop("SequenceType", None)

    decoy_str = single.pop("Decoy", None)
    decoy = decoy_str.lower() in ("true", "1", "yes") if decoy_str is not None else None

    has_annotation_identifiers = single.pop("HasAnnotationIdentifiers", "false").lower() in ("true", "1", "yes")
    proteoform_db = single.pop("ProteoformDb", "false").lower() in ("true", "1", "yes")

    return DatabaseHeader(
        prefix=prefix,
        db_name=db_name,
        db_description=db_description,
        db_version=db_version,
        db_date=db_date,
        db_sources=db_sources,
        number_of_entries=number_of_entries,
        sequence_type=sequence_type,
        decoy=decoy,
        conversion=conversion,
        has_annotation_identifiers=has_annotation_identifiers,
        proteoform_db=proteoform_db,
        custom_key_defs=custom_key_defs,
        optional_tag_defs=optional_tag_defs,
        extra=single,
    )


# ---------------------------------------------------------------------------
# Entry parsing
# ---------------------------------------------------------------------------

_DESC_PREFIX_RE = re.compile(r"^>(\S+?):(\S+)\s*(.*)")


def _parse_entry(
    description: str,
    seq_lines: list[str],
    *,
    line_no: int | None = None,
    custom_key_defs: Mapping[str, CustomKeyDef] | None = None,
) -> SequenceEntry:
    """Parse a description line and sequence lines into a SequenceEntry."""
    m = _DESC_PREFIX_RE.match(description)
    if not m:
        raise PeffParseError("Invalid description line format", line=line_no, context=description)

    prefix = m.group(1)
    db_unique_id = m.group(2)
    rest = m.group(3)

    raw_keys = split_description_keys(rest)
    sequence = "".join(sl.strip() for sl in seq_lines)

    # Defaults
    entry_id: str | None = None
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
    custom_values: dict[str, tuple[CustomKeyValue, ...]] = {}
    extra: dict[str, str] = {}

    for key, value in raw_keys.items():
        match key:
            case "ID":
                entry_id = value
            case "DbUniqueId":
                db_unique_id_key = value
            case "PName":
                pname = value
            case "GName":
                gname = value
            case "NcbiTaxId" | "OX":
                try:
                    ncbi_tax_id = int(value)
                except ValueError as err:
                    raise PeffParseError(
                        f"Invalid integer for {key}: {value!r}",
                        line=line_no,
                        context=value,
                    ) from err
            case "TaxName":
                tax_name = value
            case "Length":
                try:
                    length = int(value)
                except ValueError as err:
                    raise PeffParseError(
                        f"Invalid integer for Length: {value!r}",
                        line=line_no,
                        context=value,
                    ) from err
            case "SV":
                try:
                    sv = int(value)
                except ValueError as err:
                    raise PeffParseError(
                        f"Invalid integer for SV: {value!r}",
                        line=line_no,
                        context=value,
                    ) from err
            case "EV":
                try:
                    ev = int(value)
                except ValueError as err:
                    raise PeffParseError(
                        f"Invalid integer for EV: {value!r}",
                        line=line_no,
                        context=value,
                    ) from err
            case "PE":
                try:
                    pe = int(value)
                except ValueError as err:
                    raise PeffParseError(
                        f"Invalid integer for PE: {value!r}",
                        line=line_no,
                        context=value,
                    ) from err
            case "Decoy":
                decoy = value.lower() in ("true", "1", "yes")
            case "Comment":
                comment = value
            case "VariantSimple":
                variant_simple = _parse_variant_simple(value)
            case "VariantComplex":
                variant_complex = _parse_variant_complex(value)
            case "Variant":
                warnings.warn(
                    r"\Variant= is deprecated since PEFF 2015; use \VariantSimple= or \VariantComplex=",
                    DeprecationWarning,
                    stacklevel=2,
                )
                extra[key] = value
            case "ModResUnimod":
                mod_res_unimod = _parse_mod_res_unimod(value)
            case "ModResPsi":
                mod_res_psi = _parse_mod_res_psi(value)
            case "ModRes":
                mod_res = _parse_mod_res(value)
            case "Processed":
                processed = _parse_processed(value)
            case "DisulfideBond":
                disulfide_bond = _parse_disulfide_bond(value)
            case "Proteoform":
                proteoform = _parse_proteoform(value)
            case _:
                if custom_key_defs is not None and key in custom_key_defs:
                    custom_values[key] = _parse_custom_value(value, custom_key_defs[key])
                else:
                    extra[key] = value

    if length is not None and len(sequence) != length:
        warnings.warn(
            f"Entry {db_unique_id!r}: Length={length} but sequence has {len(sequence)} residues",
            UserWarning,
            stacklevel=2,
        )

    return SequenceEntry(
        prefix=prefix,
        db_unique_id=db_unique_id,
        sequence=sequence,
        id=entry_id,
        db_unique_id_key=db_unique_id_key,
        pname=pname,
        gname=gname,
        ncbi_tax_id=ncbi_tax_id,
        tax_name=tax_name,
        length=length,
        sv=sv,
        ev=ev,
        pe=pe,
        decoy=decoy,
        comment=comment,
        variant_simple=variant_simple,
        variant_complex=variant_complex,
        mod_res_unimod=mod_res_unimod,
        mod_res_psi=mod_res_psi,
        mod_res=mod_res,
        processed=processed,
        disulfide_bond=disulfide_bond,
        proteoform=proteoform,
        custom_values=custom_values,
        extra=extra,
    )


# ---------------------------------------------------------------------------
# PeffReader
# ---------------------------------------------------------------------------


class PeffReader:
    """Lazy reader for PEFF files."""

    def __init__(self, source: str | Path | IO[str]) -> None:
        if isinstance(source, (str, Path)):
            path = Path(source)
            self._owned_file: IO[str] | None = path.open()
            self._lines: Iterator[str] = iter(self._owned_file)
        else:
            self._owned_file = None
            self._lines = iter(source)

        self._header: FileHeader | None = None
        self._remaining: Iterator[str] | None = None
        self._defs_by_prefix: dict[str, dict[str, CustomKeyDef]] = {}

    def _ensure_header(self) -> None:
        if self._header is None:
            self._header, self._remaining = _parse_file_header(self._lines)
            for db in self._header.databases:
                if db.prefix and db.custom_key_defs:
                    self._defs_by_prefix[db.prefix] = {ckd.key_name: ckd for ckd in db.custom_key_defs}

    @property
    def header(self) -> FileHeader:
        """Parse and return the file header (cached after first access)."""
        self._ensure_header()
        assert self._header is not None
        return self._header

    def __iter__(self) -> Iterator[SequenceEntry]:
        """Yield sequence entries after the header."""
        self._ensure_header()
        assert self._remaining is not None

        current_desc: str | None = None
        current_desc_line: int | None = None
        seq_lines: list[str] = []
        line_no = 0

        for raw_line in self._remaining:
            line_no += 1
            line = raw_line.rstrip("\n\r")

            if line.startswith(">"):
                # Emit previous entry if any
                if current_desc is not None:
                    yield _parse_entry(
                        current_desc,
                        seq_lines,
                        line_no=current_desc_line,
                        custom_key_defs=self._defs_for(current_desc),
                    )
                current_desc = line
                current_desc_line = line_no
                seq_lines = []
            elif line.strip():
                seq_lines.append(line)

        # Emit last entry
        if current_desc is not None:
            yield _parse_entry(
                current_desc,
                seq_lines,
                line_no=current_desc_line,
                custom_key_defs=self._defs_for(current_desc),
            )

    def _defs_for(self, description: str) -> dict[str, CustomKeyDef] | None:
        """Look up the custom-key def map for the database matching the entry prefix."""
        if not self._defs_by_prefix:
            return None
        m = _DESC_PREFIX_RE.match(description)
        if m is None:
            return None
        return self._defs_by_prefix.get(m.group(1))

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: object) -> None:
        if self._owned_file is not None:
            self._owned_file.close()


def read_peff(source: str | Path | IO[str]) -> tuple[FileHeader, list[SequenceEntry]]:
    """Convenience: parse an entire PEFF file into header + list of entries."""
    with PeffReader(source) as reader:
        return reader.header, list(reader)
