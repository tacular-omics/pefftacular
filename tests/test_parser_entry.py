"""Tests for PEFF entry parsing."""

from io import StringIO
from pathlib import Path

from pefftacular._parser import PeffReader

FIXTURES = Path(__file__).parent / "fixtures"


class TestMinimalEntries:
    def test_entry_count(self):
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            entries = list(reader)
        assert len(entries) == 2

    def test_first_entry(self):
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            entries = list(reader)
        e = entries[0]
        assert e.prefix == "sp"
        assert e.db_unique_id == "P12345"
        assert e.pname == "Albumin"
        assert e.gname == "ALB"
        assert e.length == 10
        assert e.sequence == "MKTLLLTLVV"

    def test_second_entry_minimal(self):
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            entries = list(reader)
        e = entries[1]
        assert e.pname is None
        assert e.gname is None
        assert e.length == 5
        assert e.sequence == "ACDEF"


class TestComplexEntry:
    def _get_entry(self):
        with PeffReader(FIXTURES / "complex.peff") as reader:
            return list(reader)[0]

    def test_basic_fields(self):
        e = self._get_entry()
        assert e.prefix == "cx"
        assert e.db_unique_id == "ENTRY1"
        assert e.pname == "Test protein kinase"
        assert e.gname == "TPK1"
        assert e.ncbi_tax_id == 9606
        assert e.tax_name == "Homo sapiens"
        assert e.length == 20
        assert e.sv == 3
        assert e.ev == 1
        assert e.pe == 1

    def test_variant_simple(self):
        e = self._get_entry()
        assert len(e.variant_simple) == 2
        assert e.variant_simple[0].position == 5
        assert e.variant_simple[0].new_amino_acid == "A"
        assert e.variant_simple[0].tag is None
        assert e.variant_simple[1].position == 10
        assert e.variant_simple[1].new_amino_acid == "G"
        assert e.variant_simple[1].tag == "rs12345"

    def test_variant_complex(self):
        e = self._get_entry()
        assert len(e.variant_complex) == 1
        v = e.variant_complex[0]
        assert v.start_pos == 3
        assert v.end_pos == 7
        assert v.new_sequence == "AAALLL"
        assert v.tag == "complex_var"

    def test_mod_res_unimod(self):
        e = self._get_entry()
        assert len(e.mod_res_unimod) == 1
        m = e.mod_res_unimod[0]
        assert m.positions == (4, 7)
        assert m.accession == "UNIMOD:21"
        assert m.name == "Phospho"

    def test_mod_res_psi(self):
        e = self._get_entry()
        assert len(e.mod_res_psi) == 2
        assert e.mod_res_psi[0].positions == (6,)
        assert e.mod_res_psi[0].accession == "MOD:00048"
        assert e.mod_res_psi[1].positions == (8,)

    def test_mod_res_nested_parens(self):
        e = self._get_entry()
        assert len(e.mod_res) == 1
        m = e.mod_res[0]
        assert m.positions == (3,)
        assert m.accession == ""
        assert "GlcNAc" in m.name

    def test_processed(self):
        e = self._get_entry()
        assert len(e.processed) == 1
        p = e.processed[0]
        assert p.start_pos == 1
        assert p.end_pos == 20
        assert p.accession == "PEFF:0001"
        assert p.name == "mature protein"

    def test_sequence(self):
        e = self._get_entry()
        assert e.sequence == "ACDEFGHIKLMNPQRSTVWY"
        assert len(e.sequence) == 20


class TestContextManager:
    def test_reader_as_context_manager(self):
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            _ = reader.header
            entries = list(reader)
        assert len(entries) == 2

    def test_string_io(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n>t:X1 \\PName=Test \\Length=3\nABC\n"
        )
        with PeffReader(StringIO(data)) as reader:
            entries = list(reader)
        assert entries[0].sequence == "ABC"


class TestMultilineSequence:
    def test_sequence_reassembly(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
            ">t:X1 \\Length=12\nACDEFG\nHIKLMN\n"
        )
        with PeffReader(StringIO(data)) as reader:
            entries = list(reader)
        assert entries[0].sequence == "ACDEFGHIKLMN"
