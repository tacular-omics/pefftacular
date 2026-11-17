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


class TestDeprecatedVariantKey:
    def test_variant_warns_with_peff_warning(self):
        import pytest

        from pefftacular import PeffWarning, read_peff

        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
            ">t:X1 \\Variant=(1|1|K) \\Length=5\nMKTLL\n"
        )
        with pytest.warns(PeffWarning, match="deprecated"):
            _, entries = read_peff(StringIO(data))
        assert entries[0].extra["Variant"] == "(1|1|K)"


class TestPositionHelpers:
    def test_single_position_range(self):
        from pefftacular._parser import _parse_sequence_range

        r = _parse_sequence_range("7")
        assert (r.start, r.end) == (7, 7)

    def test_parse_position_alias_kept(self):
        from pefftacular._parser import _parse_position, parse_position

        assert parse_position is _parse_position
        assert parse_position("?") == "?"


class TestRepeatedDescriptionKeys:
    """Spec 3.3.3 forbids repeating a key, but real files do it.

    A key whose value is a list of items keeps every occurrence's items; any
    other repeated key keeps the last value and warns about the one it drops.
    """

    HEADER = (
        "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
        "# NumberOfEntries=1\n# SequenceType=AA\n"
    )
    CUSTOM = '# CustomKeyDef=(KeyName=Foo|Description="d"|FieldNames=f)\n'

    def _read(self, description, *, header_extra=""):
        from pefftacular import read_peff

        data = f"{self.HEADER}{header_extra}# //\n>t:X1 {description}\nMKTLL\n"
        return read_peff(StringIO(data))[1][0]

    def _read_quietly(self, description, *, header_extra=""):
        """Read and fail if anything warns."""
        import warnings

        from pefftacular import PeffWarning

        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            return self._read(description, header_extra=header_extra)

    def test_repeated_annotation_key_keeps_every_item(self):
        e = self._read_quietly(r"\VariantSimple=(1|A) \VariantSimple=(2|C)")
        assert [(v.position, v.new_amino_acid) for v in e.variant_simple] == [(1, "A"), (2, "C")]

    def test_items_keep_file_order_across_three_occurrences(self):
        e = self._read_quietly(r"\ModRes=(1||a) \ModRes=(2||b) \ModRes=(3||c)")
        assert [m.name for m in e.mod_res] == ["a", "b", "c"]

    def test_repeat_with_an_escape_uses_the_scanning_path(self):
        e = self._read_quietly(r"\ModRes=(1||a\|b) \ModRes=(2||c)")
        assert [(m.positions, m.name) for m in e.mod_res] == [((1,), "a|b"), ((2,), "c")]

    def test_repeated_custom_key_keeps_every_item(self):
        e = self._read_quietly(r"\Foo=(a) \Foo=(b)", header_extra=self.CUSTOM)
        assert [v.fields["f"] for v in e.custom_values["Foo"]] == ["a", "b"]
        # Same result as writing the two items under one key.
        one = self._read_quietly(r"\Foo=(a)(b)", header_extra=self.CUSTOM)
        assert e.custom_values == one.custom_values

    def test_repeated_scalar_key_warns_and_keeps_the_last(self):
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match=r"key 'PName' appears more than once.*'first' is dropped"):
            e = self._read(r"\PName=first \PName=second")
        assert e.pname == "second"
        assert "PName" not in e.extra

    def test_repeated_undeclared_key_warns_and_keeps_the_last(self):
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match="key 'Bar' appears more than once"):
            e = self._read(r"\Bar=(a) \Bar=(b)")
        assert e.extra["Bar"] == "(b)"

    def test_bare_custom_value_is_not_glued_together(self):
        """A value that is not a list of items cannot be concatenated."""
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match="key 'Foo' appears more than once"):
            e = self._read(r"\Foo=a \Foo=b", header_extra=self.CUSTOM)
        assert [v.fields["f"] for v in e.custom_values["Foo"]] == ["b"]

    def test_text_that_merely_ends_in_a_parenthesis_is_not_an_item_list(self):
        """1.4.0 raised a PeffParseError here; 1.3.0 read the last value (no warning)."""
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match=r"key 'Foo' appears more than once.*'see \(ref\)' is dropped"):
            e = self._read(r"\Foo=see (ref) \Foo=(two)", header_extra=self.CUSTOM)
        assert [v.fields["f"] for v in e.custom_values["Foo"]] == ["two"]

    def test_annotation_key_after_free_text_keeps_the_last_value(self):
        """The same for an annotation key: ``see (ref)(2|UNIMOD:1|Acetyl)`` cannot be read."""
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match="key 'ModResUnimod' appears more than once"):
            e = self._read(r"\ModResUnimod=see (ref) \ModResUnimod=(2|UNIMOD:1|Acetyl)")
        assert [(m.positions, m.name) for m in e.mod_res_unimod] == [((2,), "Acetyl")]

    def test_an_escaped_closing_paren_does_not_end_an_item_list(self):
        r"""``(a)\)`` ends with a literal ``)``, so it is not a list another list can follow."""
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match="key 'ModRes' appears more than once"):
            e = self._read(r"\ModRes=(1||a)\) \ModRes=(2||b)")
        assert [(m.positions, m.name) for m in e.mod_res] == [((2,), "b")]

    def test_a_real_item_list_still_merges_every_occurrence(self):
        # The point of the merge, unchanged: both occurrences are lists of items.
        e = self._read_quietly(r"\ModResUnimod=(1|UNIMOD:1|a)(2|UNIMOD:2|b) \ModResUnimod=(3|UNIMOD:3|c)")
        assert [(m.positions, m.name) for m in e.mod_res_unimod] == [((1,), "a"), ((2,), "b"), ((3,), "c")]

    def test_warning_names_the_entry(self):
        import pytest

        from pefftacular import PeffWarning

        with pytest.warns(PeffWarning, match="Entry t:X1"):
            self._read(r"\PName=first \PName=second")

    def test_psi_valid_example_keeps_both_modification_sites(self):
        """``PEFF_Tiny_Valid.peff`` repeats ``\\ModResUnimod`` (spec 3.3.3)."""
        import io
        import warnings

        from pefftacular import PeffWarning, read_peff, write_peff

        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            header, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        entry = next(e for e in entries if e.db_unique_id.startswith("gi"))
        assert [m.positions for m in entry.mod_res_unimod] == [(15,), ("?",)]

        buf = io.StringIO()
        write_peff(header, entries, buf)
        line = next(line for line in buf.getvalue().splitlines() if line.startswith(">nr:gi"))
        assert r"\ModResUnimod=(15|UNIMOD:4|Carbamidomethyl)(?|UNIMOD:4|Carbamidomethyl)" in line
