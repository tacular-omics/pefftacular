"""Tests for PEFF header parsing."""

import warnings
from io import StringIO
from pathlib import Path

import pytest

from pefftacular._parser import PeffReader
from pefftacular.errors import PeffParseError

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
