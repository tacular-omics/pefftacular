"""Roundtrip tests: parse -> write -> re-parse -> assert equal."""

from io import StringIO
from pathlib import Path

from pefftacular._parser import PeffReader, read_peff
from pefftacular._writer import write_peff

FIXTURES = Path(__file__).parent / "fixtures"


def _roundtrip(fixture_path: Path):
    """Parse a fixture, write to StringIO, re-parse, and compare."""
    header1, entries1 = read_peff(fixture_path)

    buf = StringIO()
    write_peff(header1, entries1, buf)
    buf.seek(0)

    header2, entries2 = read_peff(buf)
    return header1, entries1, header2, entries2


class TestRoundtripMinimal:
    def test_header_preserved(self):
        h1, _, h2, _ = _roundtrip(FIXTURES / "minimal.peff")
        assert h1.peff_version == h2.peff_version
        assert h1.general_comments == h2.general_comments
        assert len(h1.databases) == len(h2.databases)
        db1, db2 = h1.databases[0], h2.databases[0]
        assert db1.prefix == db2.prefix
        assert db1.db_name == db2.db_name
        assert db1.db_version == db2.db_version
        assert db1.number_of_entries == db2.number_of_entries
        assert db1.sequence_type == db2.sequence_type

    def test_entries_preserved(self):
        _, e1, _, e2 = _roundtrip(FIXTURES / "minimal.peff")
        assert len(e1) == len(e2)
        for a, b in zip(e1, e2, strict=True):
            assert a.prefix == b.prefix
            assert a.db_unique_id == b.db_unique_id
            assert a.sequence == b.sequence
            assert a.pname == b.pname
            assert a.gname == b.gname
            assert a.length == b.length


class TestRoundtripComplex:
    def test_all_annotations_preserved(self):
        _, e1, _, e2 = _roundtrip(FIXTURES / "complex.peff")
        a, b = e1[0], e2[0]
        assert a.pname == b.pname
        assert a.gname == b.gname
        assert a.ncbi_tax_id == b.ncbi_tax_id
        assert a.tax_name == b.tax_name
        assert a.sv == b.sv
        assert a.ev == b.ev
        assert a.pe == b.pe
        assert a.variant_simple == b.variant_simple
        assert a.variant_complex == b.variant_complex
        assert a.mod_res_unimod == b.mod_res_unimod
        assert a.mod_res_psi == b.mod_res_psi
        assert a.mod_res == b.mod_res
        assert a.processed == b.processed
        assert a.sequence == b.sequence


class TestRoundtripMultidb:
    def test_database_count(self):
        h1, _, h2, _ = _roundtrip(FIXTURES / "multidb.peff")
        assert len(h1.databases) == len(h2.databases)

    def test_entries_from_both_dbs(self):
        _, e1, _, e2 = _roundtrip(FIXTURES / "multidb.peff")
        assert len(e1) == len(e2) == 2
        assert e1[0].prefix == e2[0].prefix == "d1"
        assert e1[1].prefix == e2[1].prefix == "d2"


class TestReadPeffConvenience:
    def test_read_peff_returns_tuple(self):
        header, entries = read_peff(FIXTURES / "minimal.peff")
        assert header.peff_version == "1.0"
        assert len(entries) == 2
