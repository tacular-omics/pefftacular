"""Tests for PEFF header parsing."""

import warnings
from io import StringIO
from pathlib import Path

import pytest

from pefftacular._parser import PeffReader, read_peff
from pefftacular.errors import PeffParseError, PeffWarning

FIXTURES = Path(__file__).parent / "fixtures"


class TestPeffVersion:
    def test_valid_version(self):
        reader = PeffReader(FIXTURES / "minimal.peff")
        assert reader.header.peff_version == "1.0"

    def test_missing_version_line(self):
        with pytest.raises(PeffParseError, match="First line"):
            PeffReader(StringIO("not a header\n")).header  # noqa: B018

    def test_empty_file(self):
        with pytest.raises(PeffParseError, match="Empty file"):
            PeffReader(StringIO("")).header  # noqa: B018

    def test_unsupported_version_warns(self):
        data = (
            "# PEFF 2.0\n# //\n# Prefix=x\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=0\n# SequenceType=AA\n# //\n"
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reader = PeffReader(StringIO(data))
            _ = reader.header
            assert any("2.0" in str(warning.message) for warning in w)


class TestGeneralComments:
    def test_single_comment(self):
        reader = PeffReader(FIXTURES / "complex.peff")
        assert reader.header.general_comments == ("Complex test file with all annotation types",)

    def test_no_comments(self):
        reader = PeffReader(FIXTURES / "minimal.peff")
        assert reader.header.general_comments == ()


class TestDatabaseHeaders:
    def test_single_database(self):
        reader = PeffReader(FIXTURES / "minimal.peff")
        assert len(reader.header.databases) == 1
        db = reader.header.databases[0]
        assert db.prefix == "sp"
        assert db.db_name == "testdb"
        assert db.db_version == "2024-01"
        assert db.db_sources == ("http://example.com",)
        assert db.number_of_entries == 2
        assert db.sequence_type == "AA"

    def test_multiple_databases(self):
        reader = PeffReader(FIXTURES / "multidb.peff")
        assert len(reader.header.databases) == 2
        assert reader.header.databases[0].prefix == "d1"
        assert reader.header.databases[0].db_name == "database-one"
        assert reader.header.databases[1].prefix == "d2"
        assert reader.header.databases[1].db_name == "database-two"

    def test_missing_mandatory_key_warns(self):
        data = "# PEFF 1.0\n# //\n# Prefix=x\n# //\n"
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            reader = PeffReader(StringIO(data))
            _ = reader.header
            warning_messages = [str(warning.message) for warning in w]
            assert any("DbVersion" in msg for msg in warning_messages)
            assert any("NumberOfEntries" in msg for msg in warning_messages)


class TestHeaderCaching:
    def test_header_cached(self):
        reader = PeffReader(FIXTURES / "minimal.peff")
        h1 = reader.header
        h2 = reader.header
        assert h1 is h2


_DB_LINES = (
    "# DbName=db\n# Prefix=my\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n"
    "# CustomKeyDef=(KeyName=Foo|FieldNames=value|FieldTypes=xsd:string)\n"
)
_ENTRY = ">my:P1 \\Foo=(bar)\nPEPTIDE\n"


def _read_recording(data: str):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        header, entries = read_peff(StringIO(data))
    return header, entries, [str(x.message) for x in w if issubclass(x.category, PeffWarning)]


class TestHeaderBlockBoundaries:
    def test_missing_closing_separator_keeps_last_block(self):
        header, entries, msgs = _read_recording("# PEFF 1.0\n# //\n" + _DB_LINES + _ENTRY)
        assert [db.prefix for db in header.databases] == ["my"]
        assert header.databases[0].custom_key_defs[0].key_name == "Foo"
        assert "Foo" in entries[0].custom_values
        assert "Foo" not in entries[0].extra
        assert any("# //" in m for m in msgs)

    def test_separator_with_trailing_whitespace(self):
        data = "# PEFF 1.0\n# //  \n" + _DB_LINES + "# // \t\n" + _ENTRY
        header, entries, msgs = _read_recording(data)
        assert [db.prefix for db in header.databases] == ["my"]
        assert "Foo" in entries[0].custom_values
        assert not msgs

    def test_blank_line_inside_header_is_skipped_with_warning(self):
        data = "# PEFF 1.0\n# //\n\n" + _DB_LINES + "   \n# //\n" + _ENTRY
        header, entries, msgs = _read_recording(data)
        assert [db.prefix for db in header.databases] == ["my"]
        assert "Foo" in entries[0].custom_values
        assert sum("Blank line" in m for m in msgs) == 2  # one per blank run

    def test_blank_lines_between_header_and_entries_do_not_warn(self):
        data = "# PEFF 1.0\n# //\n" + _DB_LINES + "# //\n\n\n" + _ENTRY
        header, entries, msgs = _read_recording(data)
        assert [db.prefix for db in header.databases] == ["my"]
        assert "Foo" in entries[0].custom_values
        assert not msgs

    def test_line_numbers_after_blank_lines(self):
        data = "# PEFF 1.0\n# //\n" + _DB_LINES + "# //\n\n>my:P1 \\Length=abc\nPEPTIDE\n"
        with pytest.raises(PeffParseError) as exc:
            read_peff(StringIO(data))
        assert exc.value.line == 12
