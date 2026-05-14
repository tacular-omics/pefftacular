"""Tests for header-declared custom keys and their structural parsing."""

from __future__ import annotations

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
    SequenceEntry,
    read_peff,
    write_peff,
)
from pefftacular._lexer import split_fields

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Lexer: quote-aware split_fields
# ---------------------------------------------------------------------------


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
        reader = PeffReader(FIXTURES / "custom_keys.peff")
        defs = reader.header.databases[0].custom_key_defs
        assert len(defs) == 2
        names = [d.key_name for d in defs]
        assert names == ["SecondaryStructure", "Score"]

    def test_concept_curie_captured(self):
        reader = PeffReader(FIXTURES / "custom_keys.peff")
        ss = reader.header.databases[0].custom_key_defs[0]
        assert ss.concept_curie == "BAO:0000014"
        assert ss.field_names == ("StartPosition", "EndPosition", "CURIE", "Description")
        assert ss.field_types == ("integer", "integer", "string", "string")

    def test_regexp_with_pipes_round_trips(self):
        reader = PeffReader(FIXTURES / "custom_keys.peff")
        ss = reader.header.databases[0].custom_key_defs[0]
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
            "(KeyName=Foo|Description=\"d\"|FieldNames=d,t|FieldTypes=date,time)",
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
                "# CustomKeyDef=(KeyName=Foo|Description=\"e\"|"
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
                "(KeyName=Foo|Description=\"x\"|FieldNames=n|FieldTypes=integer)",
                "(not_a_number)",
            )
        assert kv.fields["n"] == "not_a_number"
        assert any("cannot coerce" in str(w.message) for w in caught)

    def test_regex_mismatch_warns(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            kv = self._read(
                "(KeyName=Foo|Description=\"x\"|RegExp=\"([0-9]+)\"|FieldNames=n|FieldTypes=integer)",
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
