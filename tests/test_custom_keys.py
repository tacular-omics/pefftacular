"""Tests for header-declared custom keys and their structural parsing."""

from __future__ import annotations

import dataclasses
import warnings
from io import StringIO
from pathlib import Path

import pytest

from pefftacular import (
    CustomKeyDef,
    CustomKeyValue,
    DatabaseHeader,
    FileHeader,
    PeffReader,
    PeffWarning,
    SequenceEntry,
    read_peff,
    write_peff,
)
from pefftacular._lexer import split_fields

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Lexer: quote-aware split_fields
# ---------------------------------------------------------------------------


def _header(source):
    """Parse only the header of *source* inside the reader's context manager."""
    with PeffReader(source) as reader:
        return reader.header


class TestSplitFieldsQuoteAware:
    def test_pipes_inside_quotes_are_preserved(self):
        # A RegExp value containing escaped pipes at depth 0 must not be split.
        item = r'KeyName=Foo|RegExp="([0-9]+)\|([0-9]+)\|.*"|FieldNames=a,b'
        fields = split_fields(item)
        assert fields == [
            "KeyName=Foo",
            r'RegExp="([0-9]+)\|([0-9]+)\|.*"',
            "FieldNames=a,b",
        ]

    def test_escaped_quote_inside_quoted_value(self):
        item = r'Description="quoted \"thing\" here"|other=1'
        fields = split_fields(item)
        assert fields == [r'Description="quoted \"thing\" here"', "other=1"]


# ---------------------------------------------------------------------------
# Header: CustomKeyDef parsing
# ---------------------------------------------------------------------------


class TestCustomKeyDefHeader:
    def test_multiple_custom_key_defs_preserved(self):
        header = _header(FIXTURES / "custom_keys.peff")
        defs = header.databases[0].custom_key_defs
        assert len(defs) == 2
        names = [d.key_name for d in defs]
        assert names == ["SecondaryStructure", "Score"]

    def test_concept_curie_captured(self):
        header = _header(FIXTURES / "custom_keys.peff")
        ss = header.databases[0].custom_key_defs[0]
        assert ss.concept_curie == "BAO:0000014"
        assert ss.field_names == ("StartPosition", "EndPosition", "CURIE", "Description")
        assert ss.field_types == ("integer", "integer", "string", "string")

    def test_regexp_with_pipes_round_trips(self):
        header = _header(FIXTURES / "custom_keys.peff")
        ss = header.databases[0].custom_key_defs[0]
        assert ss.regexp == r"([0-9]+)\|([0-9]+)\|([A-Za-z]+:[A-Za-z0-9]+)?\|(.+)"


# ---------------------------------------------------------------------------
# Entry: structural custom-value parsing
# ---------------------------------------------------------------------------


class TestCustomValueParsing:
    def test_regexp_driven_parse_with_coercion(self):
        _, entries = read_peff(FIXTURES / "custom_keys.peff")
        ss = entries[0].custom_values["SecondaryStructure"]
        assert len(ss) == 2
        first = ss[0]
        assert isinstance(first, CustomKeyValue)
        assert first.fields["StartPosition"] == 10
        assert first.fields["EndPosition"] == 20
        assert first.fields["CURIE"] == "ncithesaurus:C47937"
        assert first.fields["Description"] == "Helix"
        # Second item: empty CURIE group is omitted, not coerced.
        second = ss[1]
        assert second.fields["StartPosition"] == 25
        assert second.fields["EndPosition"] == 40
        assert "CURIE" not in second.fields
        assert second.fields["Description"] == "Sheet"

    def test_pipe_split_fallback_with_types(self):
        _, entries = read_peff(FIXTURES / "custom_keys.peff")
        score = entries[0].custom_values["Score"]
        assert len(score) == 1
        f = score[0].fields
        assert f["method"] == "blast"
        assert f["value"] == pytest.approx(0.99)
        assert f["passed"] is True

    def test_entry_without_custom_keys_has_empty_dict(self):
        _, entries = read_peff(FIXTURES / "custom_keys.peff")
        assert entries[1].custom_values == {}

    def test_declared_key_does_not_appear_in_extra(self):
        _, entries = read_peff(FIXTURES / "custom_keys.peff")
        assert "SecondaryStructure" not in entries[0].extra
        assert "Score" not in entries[0].extra

    def test_undeclared_key_still_goes_to_extra(self):
        data = (
            "# PEFF 1.0\n"
            "# //\n"
            "# DbName=x\n# Prefix=sp\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
            ">sp:P1 \\NotDeclared=hello\n"
            "AC\n"
        )
        _, entries = read_peff(StringIO(data))
        assert entries[0].extra == {"NotDeclared": "hello"}


# ---------------------------------------------------------------------------
# Coercion edge cases
# ---------------------------------------------------------------------------


class TestCustomValueCoercion:
    def _read(self, header_extra: str, value: str) -> CustomKeyValue:
        data = (
            "# PEFF 1.0\n"
            "# //\n# DbName=x\n# Prefix=sp\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n"
            f"# CustomKeyDef={header_extra}\n# //\n"
            f">sp:P1 \\Foo={value}\n"
            "AC\n"
        )
        _, entries = read_peff(StringIO(data))
        return entries[0].custom_values["Foo"][0]

    def test_date_and_time(self):
        from datetime import date, time

        kv = self._read(
            '(KeyName=Foo|Description="d"|FieldNames=d,t|FieldTypes=date,time)',
            "(2024-05-13|09:30:00)",
        )
        assert kv.fields["d"] == date(2024, 5, 13)
        assert kv.fields["t"] == time(9, 30, 0)

    def test_enumeration_membership_warns_out_of_set(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            data = (
                "# PEFF 1.0\n# //\n# DbName=x\n# Prefix=sp\n# DbVersion=1\n"
                "# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n"
                '# CustomKeyDef=(KeyName=Foo|Description="e"|'
                "FieldNames=v|FieldTypes=enumeration(a|b|c))\n# //\n"
                ">sp:P1 \\Foo=(z)\nAC\n"
            )
            _, entries = read_peff(StringIO(data))
        msgs = [str(w.message) for w in caught]
        assert any("enumeration" in m for m in msgs)
        # Out-of-set values are still preserved verbatim.
        assert entries[0].custom_values["Foo"][0].fields["v"] == "z"

    def test_coercion_failure_warns_and_keeps_string(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kv = self._read(
                '(KeyName=Foo|Description="x"|FieldNames=n|FieldTypes=integer)',
                "(not_a_number)",
            )
        assert kv.fields["n"] == "not_a_number"
        assert any("cannot coerce" in str(w.message) for w in caught)

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            ("true", True),
            ("True", True),
            ("TRUE", True),
            ("1", True),
            ("yes", True),
            ("YES", True),
            ("false", False),
            ("False", False),
            ("FALSE", False),
            ("0", False),
            ("no", False),
            ("NO", False),
        ],
    )
    def test_boolean_recognised_spellings(self, value, expected):
        kv = self._read(
            '(KeyName=Foo|Description="b"|FieldNames=v|FieldTypes=boolean)',
            f"({value})",
        )
        assert kv.fields["v"] is expected

    def test_boolean_unrecognised_warns_and_keeps_string(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kv = self._read(
                '(KeyName=Foo|Description="b"|FieldNames=v|FieldTypes=boolean)',
                "(maybe)",
            )
        assert kv.fields["v"] == "maybe"
        assert any("cannot coerce" in str(w.message) for w in caught)

    def test_regex_mismatch_warns(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kv = self._read(
                '(KeyName=Foo|Description="x"|RegExp="([0-9]+)"|FieldNames=n|FieldTypes=integer)',
                "(abc)",
            )
        assert kv.fields == {}
        assert kv.raw == "abc"
        assert any("does not match RegExp" in str(w.message) for w in caught)


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------


class TestCustomKeysRoundTrip:
    def test_fixture_round_trip_preserves_custom_values(self):
        header, entries = read_peff(FIXTURES / "custom_keys.peff")
        buf = StringIO()
        write_peff(header, entries, buf)
        out = buf.getvalue()
        # Header survives.
        assert "# CustomKeyDef=" in out
        assert "ConceptCURIE=BAO:0000014" in out
        # Entry annotations survive verbatim via the raw cache.
        assert r"\SecondaryStructure=(10|20|ncithesaurus:C47937|Helix)(25|40||Sheet)" in out
        assert r"\Score=(blast|0.99|true)" in out

        # Reparse and check structured fields still match.
        header2, entries2 = read_peff(StringIO(out))
        assert entries2[0].custom_values["SecondaryStructure"][0].fields["StartPosition"] == 10
        assert entries2[0].custom_values["Score"][0].fields["passed"] is True

    def test_edited_fields_are_written_not_stale_raw(self):
        header, entries = read_peff(FIXTURES / "custom_keys.peff")
        old = entries[0].custom_values["Score"][0]
        edited = dataclasses.replace(old, fields={**old.fields, "passed": False})
        entry = dataclasses.replace(entries[0], custom_values={**entries[0].custom_values, "Score": (edited,)})
        buf = StringIO()
        write_peff(header, [entry, *entries[1:]], buf)
        out = buf.getvalue()
        assert r"\Score=(blast|0.99|false)" in out
        _, entries2 = read_peff(StringIO(out))
        assert entries2[0].custom_values["Score"][0].fields["passed"] is False
        # Unedited values keep their original text.
        assert r"\SecondaryStructure=(10|20|ncithesaurus:C47937|Helix)(25|40||Sheet)" in out

    def test_unparseable_raw_is_rebuilt_from_fields(self):
        header, entries = read_peff(FIXTURES / "custom_keys.peff")
        bad = CustomKeyValue(key_name="Score", fields={"method": "m", "value": 1.5, "passed": True}, raw="a(b")
        entry = dataclasses.replace(entries[0], custom_values={"Score": (bad,)})
        buf = StringIO()
        write_peff(header, [entry, *entries[1:]], buf)
        assert r"\Score=(m|1.5|true)" in buf.getvalue()

    def test_raw_for_undeclared_key_is_written_verbatim(self):
        header, entries = read_peff(FIXTURES / "custom_keys.peff")
        v = CustomKeyValue(key_name="Undeclared", fields={"x": "other"}, raw="kept|as|is")
        entry = dataclasses.replace(entries[0], custom_values={"Undeclared": (v,)})
        buf = StringIO()
        write_peff(header, [entry, *entries[1:]], buf)
        assert r"\Undeclared=(kept|as|is)" in buf.getvalue()

    def test_write_from_constructed_fields_when_raw_empty(self):
        header = FileHeader(
            peff_version="1.0",
            databases=(
                DatabaseHeader(
                    prefix="sp",
                    db_name="x",
                    db_version="1",
                    db_sources=("x",),
                    number_of_entries=1,
                    sequence_type="AA",
                    custom_key_defs=(
                        CustomKeyDef(
                            key_name="Score",
                            description="d",
                            field_names=("method", "value"),
                            field_types=("string", "decimal"),
                        ),
                    ),
                ),
            ),
        )
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="P1",
            sequence="AC",
            custom_values={
                "Score": (
                    CustomKeyValue(
                        key_name="Score",
                        fields={"method": "blast", "value": 0.5},
                        raw="",
                    ),
                ),
            },
        )
        buf = StringIO()
        write_peff(header, [entry], buf)
        assert r"\Score=(blast|0.5)" in buf.getvalue()

    def test_constructed_fields_with_separators_roundtrip(self):
        header = FileHeader(
            peff_version="1.0",
            databases=(
                DatabaseHeader(
                    prefix="sp",
                    db_name="x",
                    db_version="1",
                    db_sources=("x",),
                    number_of_entries=1,
                    sequence_type="AA",
                    custom_key_defs=(
                        CustomKeyDef(key_name="Note", description="d", field_names=("A", "B"), field_types=()),
                    ),
                ),
            ),
        )
        hostile = r"y|z \(w) v) \ID=evil"
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="P1",
            sequence="MK",
            custom_values={"Note": (CustomKeyValue(key_name="Note", fields={"A": "x", "B": hostile}, raw=""),)},
        )
        buf = StringIO()
        write_peff(header, [entry], buf)
        buf.seek(0)
        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            _, entries = read_peff(buf)
        assert entries[0].custom_values["Note"][0].fields == {"A": "x", "B": hostile}
        assert entries[0].id is None


# ---------------------------------------------------------------------------
# A quoted CustomKeyDef value keeps its backslashes (spec 3.3.2)
# ---------------------------------------------------------------------------

_RX_HEADER = (
    "# PEFF 1.0\n"
    "# //\n"
    "# DbName=test\n"
    "# Prefix=cu\n"
    "# DbVersion=1\n"
    "# DbSource=regexp-escaping\n"
    "# SequenceType=AA\n"
    "# NumberOfEntries=1\n"
)
# The spec writes one backslash before each pipe inside the quoted value.
_RX_KEYDEF = '# CustomKeyDef=(KeyName=Score|Description="a score"|RegExp="([a-z]+\\|[0-9.]+)")\n'
_RX_ENTRY = ">cu:P1 \\Score=(blast|0.99)\nACDEF\n"


def _keydef_line(text: str) -> str:
    return next(line for line in text.splitlines() if line.startswith("# CustomKeyDef="))


def _written(text: str) -> str:
    header, entries = read_peff(StringIO(text))
    buf = StringIO()
    write_peff(header, entries, buf)
    return buf.getvalue()


class TestQuotedKeyDefBackslashes:
    def test_spec_form_round_trips_byte_identical(self) -> None:
        text = _RX_HEADER + _RX_KEYDEF + "# //\n" + _RX_ENTRY
        assert _keydef_line(_written(text)) + "\n" == _RX_KEYDEF

    def test_the_written_regexp_still_means_what_it_said(self) -> None:
        # The point of the single backslash: another tool reads the quoted value literally.
        # Doubled, `([a-z]+\\|[0-9.]+)` is an alternation of `[a-z]+\` and `[0-9.]+`, and
        # matches neither the whole value nor the pipe.
        import re

        text = _RX_HEADER + _RX_KEYDEF + "# //\n" + _RX_ENTRY
        line = _keydef_line(_written(text))
        regexp = line.split('RegExp="', 1)[1].rsplit('"', 1)[0]
        assert re.compile(regexp).fullmatch("blast|0.99") is not None

    def test_a_doubled_value_is_rewritten_in_spec_form(self) -> None:
        # What pefftacular wrote up to 1.3.0. It still reads back as the same regex, and is
        # now written the way the spec writes it.
        doubled = _RX_KEYDEF.replace("\\|", "\\\\|")
        assert doubled != _RX_KEYDEF
        text = _RX_HEADER + doubled + "# //\n" + _RX_ENTRY
        written = _written(text)
        assert _keydef_line(written) + "\n" == _RX_KEYDEF
        assert read_peff(StringIO(written))[0].databases[0].custom_key_defs[0].regexp == r"([a-z]+\|[0-9.]+)"

    @pytest.mark.parametrize(
        "value",
        [r"([a-z]+\|[0-9.]+)", "a\\\\b", "a\\", 'a"b', r"\d+", "plain"],
        ids=["spec-pipe", "two-backslashes", "trailing-backslash", "quote", "digit-class", "plain"],
    )
    def test_every_backslash_shape_survives_a_round_trip(self, value: str) -> None:
        # Values the quoted form cannot hold literally fall back to escaping every
        # backslash, so none is refused and none is written as something else.
        text = _RX_HEADER + _RX_KEYDEF + "# //\n>cu:P1\nACDEF\n"
        header, entries = read_peff(StringIO(text))
        db = header.databases[0]
        ckd = dataclasses.replace(db.custom_key_defs[0], regexp=value, description=value)
        edited = dataclasses.replace(header, databases=(dataclasses.replace(db, custom_key_defs=(ckd,)),))
        buf = StringIO()
        write_peff(edited, entries, buf)
        buf.seek(0)
        back = read_peff(buf)[0].databases[0].custom_key_defs[0]
        assert (back.regexp, back.description) == (value, value)
