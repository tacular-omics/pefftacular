"""Roundtrip tests: parse -> write -> re-parse -> assert equal."""

from io import StringIO
from pathlib import Path

import pytest

from pefftacular._parser import read_peff
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


class TestRoundtripInsulin:
    """Round-trip the insulin fixture — exercises DisulfideBond and Proteoform writers."""

    @pytest.fixture(scope="class")
    def roundtrip(self):
        return _roundtrip(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")

    def test_disulfide_bonds_preserved(self, roundtrip):
        _, e1, _, e2 = roundtrip
        a, b = e1[0], e2[0]
        assert a.disulfide_bond == b.disulfide_bond
        assert len(a.disulfide_bond) == 3
        assert a.disulfide_bond[0].annot_id_refs == (1, 2)
        assert a.disulfide_bond[0].description == "between chains"
        assert a.disulfide_bond[0].annot_id == 81

    def test_proteoforms_preserved(self, roundtrip):
        _, e1, _, e2 = roundtrip
        a, b = e1[0], e2[0]
        assert a.proteoform == b.proteoform
        assert len(a.proteoform) == 11
        assert a.proteoform[0].proteoform_id == "NX_P01308-1-pf1"
        assert a.proteoform[0].name == "preproinsulin"
        assert a.proteoform[10].annot_id_refs == (81, 82, 83)

    def test_processed_with_annot_id_preserved(self, roundtrip):
        _, e1, _, e2 = roundtrip
        a, b = e1[0], e2[0]
        assert a.processed == b.processed
        assert a.processed[0].annot_id == 77


class TestReadPeffConvenience:
    def test_read_peff_returns_tuple(self):
        header, entries = read_peff(FIXTURES / "minimal.peff")
        assert header.peff_version == "1.0"
        assert len(entries) == 2


_ESCAPE_HEADER = (
    "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
    "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
)


class TestEscapingRoundtrip:
    """Components with special characters must survive parse -> write -> re-parse (spec 3.3.3)."""

    def test_literal_pipe_in_name_roundtrips(self):
        # The escaped \| is a literal pipe within the name, not a field separator.
        data = _ESCAPE_HEADER + ">t:X1 \\ModRes=(5||name with a \\| pipe)\nACDEFACDEF\n"
        h1, e1 = read_peff(StringIO(data))
        assert e1[0].mod_res[0].name == "name with a | pipe"
        assert e1[0].mod_res[0].tag is None

        buf = StringIO()
        write_peff(h1, e1, buf)
        buf.seek(0)
        _, e2 = read_peff(buf)
        assert e1[0].mod_res == e2[0].mod_res

    def test_balanced_parens_in_name_not_escaped(self):
        # The common "N-linked (GlcNAc...)" name has balanced parens and must
        # round-trip without gaining backslashes.
        data = _ESCAPE_HEADER + ">t:X1 \\ModRes=(5||N-linked (GlcNAc...))\nACDEFACDEF\n"
        h1, e1 = read_peff(StringIO(data))
        assert e1[0].mod_res[0].name == "N-linked (GlcNAc...)"

        buf = StringIO()
        write_peff(h1, e1, buf)
        text = buf.getvalue()
        assert "N-linked (GlcNAc...)" in text  # no stray escapes
        buf.seek(0)
        _, e2 = read_peff(buf)
        assert e1[0].mod_res == e2[0].mod_res

    def test_unbalanced_paren_roundtrips(self):
        data = _ESCAPE_HEADER + ">t:X1 \\VariantSimple=(1|A|smiley :-\\) tag)\nACDEFACDEF\n"
        h1, e1 = read_peff(StringIO(data))
        assert e1[0].variant_simple[0].tag == "smiley :-) tag"
        buf = StringIO()
        write_peff(h1, e1, buf)
        buf.seek(0)
        _, e2 = read_peff(buf)
        assert e1[0].variant_simple == e2[0].variant_simple


class TestDatabaseGeneralCommentRoundtrip:
    def test_db_level_comment_preserved(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# GeneralComment=db-level note\n"
            "# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
            ">t:X1 \\Length=5\nMKTLL\n"
        )
        h1, e1 = read_peff(StringIO(data))
        assert h1.databases[0].general_comments == ("db-level note",)

        buf = StringIO()
        write_peff(h1, e1, buf)
        assert "# GeneralComment=db-level note" in buf.getvalue()
        buf.seek(0)
        h2, _ = read_peff(buf)
        assert h2.databases[0].general_comments == ("db-level note",)


class TestCaseInsensitiveFlags:
    def test_proteoformdb_uppercase_recognized(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# ProteoformDB=true\n# //\n"
            ">t:X1 \\Length=5\nMKTLL\n"
        )
        h1, _ = read_peff(StringIO(data))
        db = h1.databases[0]
        assert db.proteoform_db is True
        # Not leaked into extra under the variant spelling.
        assert "ProteoformDB" not in db.extra
