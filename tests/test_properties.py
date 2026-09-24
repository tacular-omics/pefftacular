"""Property-based tests (Hypothesis) for the PEFF reader and writer.

Properties:

1. write -> read -> write is stable: the second write gives the same text, and the
   models read back equal the models written.
2. Edited custom-key fields survive: after ``dataclasses.replace(value, fields=...)``
   on a parsed ``CustomKeyValue``, a write -> read gives the edited fields.
3. Malformed input never escapes as ``IndexError``, ``KeyError`` or any other
   non-``PeffError`` exception: the reader either returns (maybe with a
   ``PeffWarning``) or raises a ``PeffError``.

The strategies build *canonical* models, i.e. ones the format can represent exactly:

* An optional tag / description / name is ``None`` or non-empty (the writer drops an
  empty one, so ``""`` reads back as ``None``).
* Scalar free-text values (``pname``, ``gname``, ``tax_name``, ``comment``) do not end
  in whitespace: the spec separates ``\\Key=value`` pairs with spaces and has no escape
  for a space, so trailing whitespace before the next key is not recoverable.
* A ``Proteoform`` id does not look like ``<digits>:...``, which the spec reserves for
  an annotation identifier (section 3.4.2).
"""

from __future__ import annotations

import dataclasses
import string
import warnings
from datetime import date, time
from io import StringIO

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

from pefftacular import (
    CustomKeyDef,
    CustomKeyValue,
    DatabaseHeader,
    DisulfideBond,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    OptionalTagDef,
    PeffError,
    PeffWriteError,
    Processed,
    Proteoform,
    SequenceEntry,
    SequenceRange,
    VariantComplex,
    VariantSimple,
    read_peff,
    write_peff,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write(header: FileHeader, entries: list[SequenceEntry]) -> str:
    buf = StringIO()
    write_peff(header, entries, buf)
    return buf.getvalue()


def _read(text: str) -> tuple[FileHeader, list[SequenceEntry]]:
    """Read *text*, ignoring warnings: the properties are about data, not validity."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_peff(StringIO(text))


# ---------------------------------------------------------------------------
# Leaf strategies
# ---------------------------------------------------------------------------

# Any character except line breaks (and lone surrogates, which cannot be encoded).
_LINE_CHARS = st.characters(exclude_categories=("Cs",), exclude_characters="\n\r")
# Characters the lexer treats specially, to make sure escaping is exercised.
_SPECIAL = st.sampled_from(list("\\|()= ,:;?-[]\"'"))
_text_char = st.one_of(_LINE_CHARS, _SPECIAL, _SPECIAL)

component = st.text(_text_char, max_size=12)
nonempty_component = st.text(_text_char, min_size=1, max_size=12)
optional_component = st.none() | nonempty_component
scalar_text = st.text(_text_char, max_size=12).filter(lambda s: not s[-1:].isspace())

identifier = st.text(string.ascii_letters + string.digits + "_", min_size=1, max_size=8)
# Tokens without whitespace or backslash: prefixes, ids and other unescaped values.
token = st.text(string.ascii_letters + string.digits + "_-.:", min_size=1, max_size=10)
prefix = st.text(string.ascii_letters + string.digits + "_-.", min_size=1, max_size=6)
header_value = st.text(_LINE_CHARS, max_size=16)

residues = st.text("ACDEFGHIKLMNPQRSTVWYBZXJOU*-", min_size=1, max_size=150)
position = st.integers(1, 5000) | st.just("?")
annot_id = st.none() | st.integers(0, 500)
small_int = st.integers(-(10**9), 10**9)

_KNOWN_ENTRY_KEYS = {
    "ID",
    "DbUniqueId",
    "PName",
    "GName",
    "NcbiTaxId",
    "OX",
    "TaxName",
    "Length",
    "SV",
    "EV",
    "PE",
    "Decoy",
    "Comment",
    "VariantSimple",
    "VariantComplex",
    "Variant",
    "ModResUnimod",
    "ModResPsi",
    "ModRes",
    "Processed",
    "DisulfideBond",
    "Proteoform",
}

# ---------------------------------------------------------------------------
# Annotation strategies
# ---------------------------------------------------------------------------

variant_simple = st.builds(
    VariantSimple, position=position, new_amino_acid=component, tag=optional_component, annot_id=annot_id
)
variant_complex = st.builds(
    VariantComplex,
    start_pos=position,
    end_pos=position,
    new_sequence=component,
    tag=optional_component,
    annot_id=annot_id,
)
_positions = st.lists(position, min_size=1, max_size=4).map(tuple)


def _mod(cls: type) -> st.SearchStrategy:
    return st.builds(
        cls, positions=_positions, accession=component, name=component, tag=optional_component, annot_id=annot_id
    )


processed = st.builds(
    Processed,
    start_pos=position,
    end_pos=position,
    accession=component,
    name=component,
    tag=optional_component,
    annot_id=annot_id,
)
disulfide_bond = st.builds(
    DisulfideBond,
    annot_id_refs=st.lists(st.integers(0, 500), min_size=1, max_size=3).map(tuple),
    description=optional_component,
    annot_id=annot_id,
)
proteoform = st.builds(
    Proteoform,
    proteoform_id=nonempty_component.filter(lambda s: not s.split(":", 1)[0].strip().isdigit() or ":" not in s),
    ranges=st.lists(st.builds(SequenceRange, start=st.integers(1, 5000), end=st.integers(1, 5000)), max_size=3).map(
        tuple
    ),
    annot_id_refs=st.lists(st.integers(0, 500), max_size=4).map(tuple),
    name=optional_component,
    annot_id=annot_id,
)


def _tuple(strategy: st.SearchStrategy) -> st.SearchStrategy:
    return st.lists(strategy, max_size=3).map(tuple)


extra_values = st.dictionaries(
    identifier.map(lambda k: f"X{k}").filter(lambda k: k not in _KNOWN_ENTRY_KEYS),
    st.text(string.ascii_letters + string.digits + "_-.:,;", max_size=10),
    max_size=2,
)

sequence_entry = st.builds(
    SequenceEntry,
    prefix=st.just("sp"),
    db_unique_id=token,
    sequence=residues,
    id=st.none() | token,
    db_unique_id_key=st.none() | token,
    pname=st.none() | scalar_text,
    gname=st.none() | scalar_text,
    ncbi_tax_id=st.none() | small_int,
    tax_name=st.none() | scalar_text,
    length=st.none() | st.integers(0, 10_000),
    sv=st.none() | small_int,
    ev=st.none() | small_int,
    pe=st.none() | small_int,
    decoy=st.none() | st.booleans(),
    comment=st.none() | scalar_text,
    variant_simple=_tuple(variant_simple),
    variant_complex=_tuple(variant_complex),
    mod_res_unimod=_tuple(_mod(ModResUnimod)),
    mod_res_psi=_tuple(_mod(ModResPsi)),
    mod_res=_tuple(_mod(ModRes)),
    processed=_tuple(processed),
    disulfide_bond=_tuple(disulfide_bond),
    proteoform=_tuple(proteoform),
    custom_values=st.just({}),
    extra=extra_values,
)

# ---------------------------------------------------------------------------
# Header strategies
# ---------------------------------------------------------------------------

_KNOWN_HEADER_KEYS = {
    "prefix",
    "dbname",
    "dbdescription",
    "dbversion",
    "dbdate",
    "dbsource",
    "generalcomment",
    "numberofentries",
    "sequencetype",
    "decoy",
    "conversion",
    "hasannotationidentifiers",
    "hasannotationidentifier",
    "proteoformdb",
    "isproteoformdb",
    "customkeydef",
    "optionaltagdef",
    "specifickey",
    "specificvalue",
}

custom_key_def = st.builds(
    CustomKeyDef,
    key_name=identifier.map(lambda k: f"C{k}"),
    description=header_value,
    concept_curie=st.none() | token,
    regexp=st.none() | header_value,
    field_names=st.lists(identifier, max_size=3).map(tuple),
    field_types=st.lists(st.sampled_from(["string", "integer", "decimal", "boolean", "date", "time"]), max_size=3).map(
        tuple
    ),
)
optional_tag_def = st.builds(
    OptionalTagDef,
    tag=st.text(_LINE_CHARS, min_size=1, max_size=8).filter(lambda s: ":" not in s),
    description=header_value,
)


def _database_header(pfx: str) -> st.SearchStrategy:
    return st.builds(
        DatabaseHeader,
        prefix=st.just(pfx),
        db_name=st.none() | header_value,
        db_description=st.none() | header_value,
        db_version=st.none() | header_value,
        db_date=st.none() | header_value,
        db_sources=_tuple(header_value),
        general_comments=_tuple(header_value),
        number_of_entries=st.none() | st.integers(0, 10**6),
        sequence_type=st.none() | st.sampled_from(["AA", "NA"]),
        decoy=st.none() | st.booleans(),
        conversion=st.none() | header_value,
        has_annotation_identifiers=st.booleans(),
        proteoform_db=st.booleans(),
        custom_key_defs=_tuple(custom_key_def),
        optional_tag_defs=_tuple(optional_tag_def),
        extra=st.dictionaries(
            identifier.map(lambda k: f"X{k}").filter(lambda k: k.lower() not in _KNOWN_HEADER_KEYS),
            header_value,
            max_size=2,
        ),
    )


file_header = st.builds(
    FileHeader,
    peff_version=st.just("1.0"),
    general_comments=_tuple(header_value),
    databases=st.tuples(_database_header("sp"), _database_header("tr")),
)


# ---------------------------------------------------------------------------
# Property 1: write -> read -> write is stable
# ---------------------------------------------------------------------------


@given(header=file_header, entries=st.lists(sequence_entry, max_size=3))
def test_write_read_write_is_stable(header: FileHeader, entries: list[SequenceEntry]) -> None:
    text = _write(header, entries)
    header2, entries2 = _read(text)
    assert header2 == header
    assert entries2 == entries
    assert _write(header2, entries2) == text


@given(entries=st.lists(sequence_entry, min_size=1, max_size=3))
def test_read_write_is_idempotent_on_written_text(entries: list[SequenceEntry]) -> None:
    header = FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp", db_name="x"),))
    text = _write(header, entries)
    assert _write(*_read(text)) == text


# ---------------------------------------------------------------------------
# Property 2: edited custom-key fields survive a write
# ---------------------------------------------------------------------------

_FIELD_VALUES: dict[str, st.SearchStrategy] = {
    "string": component,
    "integer": st.integers(-(10**12), 10**12),
    "decimal": st.floats(allow_nan=False),
    "boolean": st.booleans(),
    "date": st.dates(),
    "time": st.times(),
    "enumeration(helix|sheet|coil)": st.sampled_from(["helix", "sheet", "coil"]),
}


@st.composite
def custom_key_case(draw: st.DrawFn) -> tuple[CustomKeyDef, dict, dict]:
    """A CustomKeyDef without a RegExp, plus an original and an edited set of fields."""
    types = draw(st.lists(st.sampled_from(sorted(_FIELD_VALUES)), min_size=1, max_size=4))
    names = tuple(f"f{i}" for i in range(len(types)))
    ckd = CustomKeyDef(key_name="Custom", description="d", field_names=names, field_types=tuple(types))
    original = {n: draw(_FIELD_VALUES[t]) for n, t in zip(names, types, strict=True)}
    edited = {n: draw(_FIELD_VALUES[t]) for n, t in zip(names, types, strict=True)}
    return ckd, original, edited


def _custom_header(ckd: CustomKeyDef) -> FileHeader:
    return FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp", custom_key_defs=(ckd,)),))


@given(case=custom_key_case())
def test_edited_custom_fields_survive(case: tuple[CustomKeyDef, dict, dict]) -> None:
    ckd, original, edited = case
    header = _custom_header(ckd)
    built = SequenceEntry(
        prefix="sp", db_unique_id="P1", sequence="MK", custom_values={"Custom": (CustomKeyValue("Custom", original),)}
    )
    # Read it once so the value carries the parsed raw text, as real input would.
    _, (parsed,) = _read(_write(header, [built]))
    (value,) = parsed.custom_values["Custom"]
    assert value.fields == original

    edited_value = dataclasses.replace(value, fields=edited)
    entry = dataclasses.replace(parsed, custom_values={"Custom": (edited_value,)})
    _, (reread,) = _read(_write(header, [entry]))
    assert reread.custom_values["Custom"][0].fields == edited


# RegExp-controlled key: StartPosition|EndPosition|Label, the shape of the spec's
# SecondaryStructure example (section 3.3.2).
_REGEXP_DEF = CustomKeyDef(
    key_name="Custom",
    description="d",
    regexp=r"([0-9]+)\|([0-9]+)\|(.*)",
    field_names=("start", "end", "label"),
    field_types=("integer", "integer", "string"),
)


@given(
    start=st.integers(0, 10**6),
    end=st.integers(0, 10**6),
    label=component,
)
def test_edited_regexp_custom_fields_survive_or_raise(start: int, end: int, label: str) -> None:
    header = _custom_header(_REGEXP_DEF)
    _, (parsed,) = _read(
        ">".join(
            [_write(header, []), "sp:P1 \\Custom=(1|2|x)\nMK\n"],
        )
    )
    (value,) = parsed.custom_values["Custom"]
    edited = {"start": start, "end": end, "label": label}
    entry = dataclasses.replace(parsed, custom_values={"Custom": (dataclasses.replace(value, fields=edited),)})
    try:
        text = _write(header, [entry])
    except PeffWriteError:
        return  # the value cannot be expressed through this RegExp: refused, not corrupted
    _, (reread,) = _read(text)
    assert reread.custom_values["Custom"][0].fields == edited


# ---------------------------------------------------------------------------
# Property 3: malformed input raises PeffError (or warns), nothing else
# ---------------------------------------------------------------------------


def _read_or_peff_error(text: str) -> None:
    try:
        _read(text)
    except PeffError:
        pass


_HEADER = (
    "# PEFF 1.0\n# //\n# DbName=x\n# Prefix=sp\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n"
    "# SequenceType=AA\n"
    '# CustomKeyDef=(KeyName=Typed|Description="d"|FieldNames=a,b,c|FieldTypes=integer,date,boolean)\n'
    '# CustomKeyDef=(KeyName=Rx|Description="d"|RegExp="([0-9]+)\\|(.*)"|FieldNames=n,s|FieldTypes=integer,time)\n'
    "# //\n"
)

_DESC_TOKENS = st.sampled_from(
    [
        *(f"\\{k}=" for k in sorted(_KNOWN_ENTRY_KEYS)),
        "\\Typed=",
        "\\Rx=",
        "\\",
        "\\\\",
        "\\|",
        "\\(",
        "\\)",
        "(",
        ")",
        ")(",
        "|",
        "||",
        ",",
        ":",
        "-",
        " ",
        "?",
        "=",
        "1",
        "42",
        "0:",
        "A",
        "MOD:00046",
        "x",
        "true",
        "2024-01-01",
        "12:30",
        "1-10",
    ]
)

_HEADER_TOKENS = st.sampled_from(
    [
        "# ",
        "#",
        "//",
        "\n",
        "DbName=",
        "Prefix=",
        "NumberOfEntries=",
        "NumberOfEntries=x",
        "CustomKeyDef=",
        "OptionalTagDef=",
        "GeneralComment=",
        "Decoy=",
        "ProteoformDB=",
        "KeyName=",
        "RegExp=",
        "FieldNames=",
        "FieldTypes=",
        "enumeration(",
        "integer",
        '"',
        '\\"',
        "(",
        ")",
        "|",
        ",",
        "=",
        "[",
        "*",
        "a",
        "1",
    ]
)


@given(st.text(max_size=300))
def test_arbitrary_text_raises_only_peff_error(text: str) -> None:
    _read_or_peff_error(text)


@given(st.lists(_DESC_TOKENS, max_size=30).map("".join), residues)
def test_malformed_description_raises_only_peff_error(desc: str, seq: str) -> None:
    _read_or_peff_error(f"{_HEADER}>sp:P1 {desc}\n{seq}\n")


@given(st.lists(_HEADER_TOKENS, max_size=40).map("".join))
def test_malformed_header_raises_only_peff_error(body: str) -> None:
    _read_or_peff_error(f"# PEFF 1.0\n# //\n# {body}\n# //\n>sp:P1 \\Typed=(1|2|3)\nMK\n")


@given(
    header=file_header,
    entries=st.lists(sequence_entry, min_size=1, max_size=2),
    data=st.data(),
)
def test_mutated_valid_file_raises_only_peff_error(
    header: FileHeader, entries: list[SequenceEntry], data: st.DataObject
) -> None:
    text = _write(header, entries)
    assume(text)
    for _ in range(data.draw(st.integers(1, 4))):
        i = data.draw(st.integers(0, len(text)))
        j = data.draw(st.integers(i, min(len(text), i + 3)))
        insert = data.draw(st.text(st.sampled_from(list("()|\\=:,>#/ \n?-1A")), max_size=3))
        text = text[:i] + insert + text[j:]
    _read_or_peff_error(text)


# ---------------------------------------------------------------------------
# Regressions found by the properties above
# ---------------------------------------------------------------------------


def test_time_and_date_fields_round_trip() -> None:
    ckd = CustomKeyDef(key_name="Custom", description="d", field_names=("t", "d"), field_types=("time", "date"))
    fields = {"t": time(1, 2, 3, 4), "d": date(9, 1, 2)}
    entry = SequenceEntry(
        prefix="sp", db_unique_id="P1", sequence="MK", custom_values={"Custom": (CustomKeyValue("Custom", fields),)}
    )
    _, (reread,) = _read(_write(_custom_header(ckd), [entry]))
    assert reread.custom_values["Custom"][0].fields == fields


def test_proteoform_id_with_paren_or_pipe_is_escaped() -> None:
    pf = (
        Proteoform(proteoform_id=")", ranges=(), annot_id_refs=()),
        Proteoform(proteoform_id="a|b", ranges=(SequenceRange(1, 2),), annot_id_refs=()),
    )
    entry = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK", proteoform=pf)
    header = FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp"),))
    _, (reread,) = _read(_write(header, [entry]))
    assert reread.proteoform == pf


def test_empty_disulfide_description_is_not_written() -> None:
    entry = SequenceEntry(
        prefix="sp", db_unique_id="P1", sequence="MK", disulfide_bond=(DisulfideBond((1, 2), description=""),)
    )
    header = FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp"),))
    text = _write(header, [entry])
    assert "\\DisulfideBond=(1,2)\n" in text
    assert _write(*_read(text)) == text


def test_edited_regexp_value_that_breaks_the_item_raises() -> None:
    header = _custom_header(_REGEXP_DEF)
    value = CustomKeyValue("Custom", {"start": 1, "end": 2, "label": ")"})
    entry = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK", custom_values={"Custom": (value,)})
    with pytest.raises(PeffWriteError, match="RegExp"):
        _write(header, [entry])


def test_custom_key_def_quoted_values_may_hold_unbalanced_parens() -> None:
    for ckd in (
        CustomKeyDef(key_name="Custom", description="left ( only"),
        CustomKeyDef(key_name="Custom", description="right ) only"),
        CustomKeyDef(key_name="Custom", description="d", regexp="([0-9]+"),
    ):
        header2, _ = _read(_write(_custom_header(ckd), []))
        assert header2.databases[0].custom_key_defs == (ckd,)


_OK = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK")
_OK_HEADER = FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp"),))


@pytest.mark.parametrize(
    ("header", "entry", "match"),
    [
        (_OK_HEADER, dataclasses.replace(_OK, pname="a\nb"), "line break"),
        (_OK_HEADER, dataclasses.replace(_OK, mod_res=(ModRes((1,), "MOD:1", "x\r"),)), "line break"),
        (_OK_HEADER, dataclasses.replace(_OK, extra={"Xk": "a\nb"}), "line break"),
        (dataclasses.replace(_OK_HEADER, general_comments=("a\nb",)), _OK, "line break"),
        (_OK_HEADER, dataclasses.replace(_OK, prefix="s p"), "prefix"),
        (_OK_HEADER, dataclasses.replace(_OK, prefix="s:p"), "prefix"),
        (_OK_HEADER, dataclasses.replace(_OK, db_unique_id="P 1"), "db_unique_id"),
        (_OK_HEADER, dataclasses.replace(_OK, sequence="MK LL"), "sequence"),
        (_OK_HEADER, dataclasses.replace(_OK, sequence="MK>LL"), "sequence"),
    ],
)
def test_values_that_would_corrupt_the_file_raise(header: FileHeader, entry: SequenceEntry, match: str) -> None:
    buf = StringIO()
    with pytest.raises(PeffWriteError, match=match):
        write_peff(header, [entry], buf)
    assert buf.getvalue() == ""
