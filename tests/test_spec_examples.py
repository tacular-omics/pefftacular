"""Every example in the PEFF 1.0 specification, parsed and written back.

Reference: ``PEFF_SpecDoc_1.0_FINAL.pdf`` (HUPO-PSI, 2019-04-05), in the repo root.
Section numbers below are the spec's. The multi-line examples are transcribed into
``tests/fixtures/spec/``: the PDF wraps long description lines, so each was rejoined
into one line (a space where the wrap is at a word break, nothing where it splits a
token such as ``(12`` / ``9|R)``); the sequences are copied unchanged. The two entry
examples are printed without a file header, so each fixture adds a minimal one.

Each legal example must parse to the documented structure without a ``PeffWarning`` and
survive write -> read unchanged. Each example the spec marks ILLEGAL must give a
``PeffWarning`` or a ``PeffParseError``, never a silent success.
"""

import warnings
from io import StringIO
from pathlib import Path

import pytest

from pefftacular import (
    CustomKeyDef,
    FileHeader,
    PeffParseError,
    PeffReader,
    PeffWarning,
    SequenceEntry,
    read_peff,
    write_peff,
)

SPEC = Path(__file__).parent / "fixtures" / "spec"

_HEADER = (
    "# PEFF 1.0\n# //\n# DbName=spec\n# Prefix=sp\n# DbVersion=1\n# DbSource=spec\n"
    "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
)
_SEQ = "ACDEFGHIKLMNPQRSTVWY" * 20  # 400 residues: every spec position fits


def _read_strict(text: str) -> tuple[FileHeader, list[SequenceEntry]]:
    """Parse *text*, turning any PeffWarning into an error."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", PeffWarning)
        return read_peff(StringIO(text))


def _read_header_strict(text: str) -> FileHeader:
    """Parse only the header of *text* (no entry-count check), PeffWarning -> error."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", PeffWarning)
        return PeffReader(StringIO(text)).header


def _entry(desc: str, seq: str = _SEQ) -> SequenceEntry:
    _, entries = _read_strict(f"{_HEADER}>sp:P1 {desc}\n{seq}\n")
    return entries[0]


def _write(header: FileHeader, entries: list[SequenceEntry]) -> str:
    buf = StringIO()
    write_peff(header, entries, buf)
    return buf.getvalue()


def _assert_round_trip(header: FileHeader, entries: list[SequenceEntry]) -> None:
    """write -> read gives equal models, and a second write gives the same text."""
    text = _write(header, entries)
    if entries:
        header2, entries2 = _read_strict(text)
    else:
        header2, entries2 = _read_header_strict(text), []
    assert header2 == header
    assert entries2 == entries
    assert _write(header2, entries2) == text


def _assert_entry_round_trip(desc: str) -> None:
    header, entries = _read_strict(f"{_HEADER}>sp:P1 {desc}\n{_SEQ}\n")
    _assert_round_trip(header, entries)


# ---------------------------------------------------------------------------
# 3.3.1 file header example
# ---------------------------------------------------------------------------


class TestHeaderExample:
    """Section 3.3.1: two database blocks, database-level GeneralComment, flags."""

    def _header(self) -> FileHeader:
        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            with PeffReader(SPEC / "spec_3_3_1_header.peff") as reader:
                return reader.header

    def test_file_description_block(self) -> None:
        header = self._header()
        assert header.peff_version == "1.0"
        assert header.general_comments == ("This is a hand-crafted example comment",)
        assert len(header.databases) == 2

    def test_nextprot_block(self) -> None:
        db = self._header().databases[0]
        assert db.db_name == "neXtProt-extract"
        assert db.prefix == "nxp"
        assert db.db_description == "extract of neXtProt with manual modifications"
        assert db.general_comments == ("A GeneralComment specific to one database is also legal here",)
        assert db.decoy is False
        assert db.db_version == "2018-01-11"
        assert db.db_sources == ("http://www.nextprot.org",)
        assert db.number_of_entries == 62
        assert db.sequence_type == "AA"
        assert db.proteoform_db is False
        assert db.has_annotation_identifiers is False
        assert db.extra == {}

    def test_mydb_block(self) -> None:
        db = self._header().databases[1]
        assert db.db_name == "myDB"
        assert db.prefix == "my"
        assert db.db_description == "FGF21 proteoforms from top-down experiment PXD123456"
        assert db.db_version == "1.1"
        assert db.db_sources == ("PXD123456",)
        assert db.number_of_entries == 2
        assert db.sequence_type == "AA"
        # The example spells it "ProteoformDB"; section 3.4.1 spells it "ProteoformDb".
        assert db.proteoform_db is True
        assert db.has_annotation_identifiers is True
        assert db.extra == {}

    def test_writes_back(self) -> None:
        _assert_round_trip(self._header(), [])


# ---------------------------------------------------------------------------
# 3.3.2 custom keys
# ---------------------------------------------------------------------------

# Rejoined from the PDF, which wraps the line after "Secondary structure",
# "(.*)" and "integer,i".
_SPEC_CUSTOM_KEY_DEF = (
    '# CustomKeyDef=(KeyName=SecondaryStructure|Description="Secondary structure information"'
    r'|ConceptCURIE=BAO:0000014|RegExp="([0-9]+)\|([0-9]+)\|([A-Za-z]+:[0-9]+)?\|(.*)\|?(.+)?"'
    "|FieldNames=StartPosition,EndPosition,CURIE,Description,OptionalTag"
    "|FieldTypes=integer,integer,string,string,string)\n"
)
_SPEC_CUSTOM_VALUE = r"\SecondaryStructure=(617|673|ncithesaurus:C47937|Helix)"


def _custom_key_file(keydef_line: str) -> str:
    return _HEADER.removesuffix("# //\n") + keydef_line + "# //\n"


class TestCustomKeyExample:
    """Section 3.3.2: the SecondaryStructure CustomKeyDef and its entry value."""

    def test_definition_parses(self) -> None:
        header = _read_header_strict(_custom_key_file(_SPEC_CUSTOM_KEY_DEF))
        assert header.databases[0].custom_key_defs == (
            CustomKeyDef(
                key_name="SecondaryStructure",
                description="Secondary structure information",
                concept_curie="BAO:0000014",
                regexp=r"([0-9]+)\|([0-9]+)\|([A-Za-z]+:[0-9]+)?\|(.*)\|?(.+)?",
                field_names=("StartPosition", "EndPosition", "CURIE", "Description", "OptionalTag"),
                field_types=("integer", "integer", "string", "string", "string"),
            ),
        )

    def test_spec_value_does_not_match_spec_regexp(self) -> None:
        # Spec erratum: the RegExp's CURIE group is [A-Za-z]+:[0-9]+, but the example
        # CURIE "ncithesaurus:C47937" has a letter after the colon, so the spec's own
        # value does not match. The reader warns and keeps the raw text.
        text = f"{_custom_key_file(_SPEC_CUSTOM_KEY_DEF)}>sp:P1 {_SPEC_CUSTOM_VALUE}\n{_SEQ}\n"
        with pytest.warns(PeffWarning, match="does not match RegExp"):
            header, entries = read_peff(StringIO(text))
        (value,) = entries[0].custom_values["SecondaryStructure"]
        assert value.fields == {}
        assert value.raw == "617|673|ncithesaurus:C47937|Helix"
        assert _SPEC_CUSTOM_VALUE in _write(header, entries)

    def test_spec_value_with_corrected_regexp(self) -> None:
        keydef = _SPEC_CUSTOM_KEY_DEF.replace("[A-Za-z]+:[0-9]+", "[A-Za-z]+:[A-Za-z0-9]+")
        header, entries = _read_strict(f"{_custom_key_file(keydef)}>sp:P1 {_SPEC_CUSTOM_VALUE}\n{_SEQ}\n")
        (value,) = entries[0].custom_values["SecondaryStructure"]
        assert value.fields == {
            "StartPosition": 617,
            "EndPosition": 673,
            "CURIE": "ncithesaurus:C47937",
            "Description": "Helix",
        }
        _assert_round_trip(header, entries)


# ---------------------------------------------------------------------------
# 3.3.3 entry section: generic rules and the real example
# ---------------------------------------------------------------------------


class TestEntryRules:
    """Section 3.3.3 rules that come with an inline example."""

    def test_embedded_parentheses(self) -> None:
        entry = _entry(r"\ModRes=(380||N-linked (GlcNAc...))")
        (mod,) = entry.mod_res
        assert (mod.positions, mod.accession, mod.name) == ((380,), "", "N-linked (GlcNAc...)")
        _assert_entry_round_trip(r"\ModRes=(380||N-linked (GlcNAc...))")

    def test_escaped_optional_tag(self) -> None:
        # 'if an optionalTag component was "Abcg2|meta\x10", it should be written as
        # "Abcg2\|meta\\x10" in PEFF.'
        desc = r"\ModResUnimod=(100|UNIMOD:21|Phospho|Abcg2\|meta\\x10)"
        entry = _entry(desc)
        assert entry.mod_res_unimod[0].tag == r"Abcg2|meta\x10"
        header, entries = _read_strict(f"{_HEADER}>sp:P1 {desc}\n{_SEQ}\n")
        assert desc in _write(header, entries)

    def test_spaces_between_items(self) -> None:
        entry = _entry(r"\ModResUnimod=(1|UNIMOD:21|Phospho) (2|UNIMOD:21|Phospho)")
        assert [m.positions for m in entry.mod_res_unimod] == [(1,), (2,)]

    def test_crlf_line_endings(self) -> None:
        # "A CR (ASCII 13) MAY precede the LF and MUST be ignored by parsers."
        text = f"{_HEADER}>sp:P1 \\PName=Some protein \\Length=4\nMKLV\n"
        _, lf = _read_strict(text)
        _, crlf = _read_strict(text.replace("\n", "\r\n"))
        assert crlf == lf
        assert crlf[0].pname == "Some protein"

    def test_blank_lines_between_entries(self) -> None:
        text = _HEADER.replace("NumberOfEntries=1", "NumberOfEntries=2") + ">sp:P1\nMK\n\n\n>sp:P2\n\nLV\n\n"
        _, entries = _read_strict(text)
        assert [(e.db_unique_id, e.sequence) for e in entries] == [("P1", "MK"), ("P2", "LV")]

    def test_splice_variant_ids(self) -> None:
        # Section 3.5.1: ">sp:P01234-1 and >sp:P01234-2".
        text = _HEADER.replace("NumberOfEntries=1", "NumberOfEntries=2") + ">sp:P01234-1\nMK\n>sp:P01234-2\nLV\n"
        _, entries = _read_strict(text)
        assert [e.db_unique_id for e in entries] == ["P01234-1", "P01234-2"]

    def test_allowed_sequence_characters(self) -> None:
        # Both residue tables, including "*" (interruption) and "-" (gap).
        seq = "ARNDCQEGHILKMFPOSUTWYVBZXJ*GATCURYKMSWBDHVN-"
        assert _entry("", seq).sequence == seq


class TestTyro3Example:
    """Section 3.3.3 "Real example": neXtProt NX_Q06418-1 (TYRO3)."""

    def _read(self) -> tuple[FileHeader, list[SequenceEntry]]:
        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            return read_peff(SPEC / "spec_3_3_3_tyro3.peff")

    def test_scalar_keys(self) -> None:
        (entry,) = self._read()[1]
        assert (entry.prefix, entry.db_unique_id) == ("nxp", "NX_Q06418-1")
        assert entry.pname == "Tyrosine-protein kinase receptor TYRO3 isoform Iso 1"
        assert entry.gname == "TYRO3"
        assert entry.ncbi_tax_id == 9606
        assert entry.tax_name == "Homo Sapiens"
        assert (entry.length, entry.sv, entry.ev, entry.pe) == (890, 135, 357, 1)
        assert len(entry.sequence) == 890
        assert entry.sequence.startswith("MALRRSMGRPGLPP")
        assert entry.sequence.endswith("LLLLQQGLLPHSSC")
        assert entry.extra == {}

    def test_processed(self) -> None:
        (entry,) = self._read()[1]
        assert [(p.start_pos, p.end_pos, p.accession, p.name) for p in entry.processed] == [
            (1, 40, "PEFF:0001021", "signal peptide"),
            (41, 890, "PEFF:0001020", "mature protein"),
        ]

    def test_mod_res_psi(self) -> None:
        (entry,) = self._read()[1]
        phospho = "O4'-phospho-L-tyrosine"
        assert [(m.positions, m.accession, m.name) for m in entry.mod_res_psi] == [
            ((681,), "MOD:00048", phospho),
            ((685,), "MOD:00048", phospho),
            ((686,), "MOD:00048", phospho),
            ((804,), "MOD:00048", phospho),
            ((64,), "MOD:00798", "half cystine"),
            ((117,), "MOD:00798", "half cystine"),
            ((160,), "MOD:00798", "half cystine"),
            ((203,), "MOD:00798", "half cystine"),
        ]

    def test_mod_res(self) -> None:
        (entry,) = self._read()[1]
        assert [m.positions[0] for m in entry.mod_res] == [63, 191, 230, 240, 293, 366, 380]
        assert {(m.accession, m.name, m.tag) for m in entry.mod_res} == {("", "N-linked (GlcNAc...)", None)}

    def test_variant_simple(self) -> None:
        (entry,) = self._read()[1]
        assert len(entry.variant_simple) == 113
        first, last = entry.variant_simple[0], entry.variant_simple[-1]
        assert (first.position, first.new_amino_acid) == (21, "L")
        assert (last.position, last.new_amino_acid) == (875, "R")
        # (12 / 9|R) is split across two PDF lines; the rejoined item is (129|R).
        assert any(v.position == 129 and v.new_amino_acid == "R" for v in entry.variant_simple)

    def test_writes_back(self) -> None:
        _assert_round_trip(*self._read())


# ---------------------------------------------------------------------------
# 3.3.8-3.3.13 annotation keys: legal examples
# ---------------------------------------------------------------------------


class TestVariantSimpleExamples:
    """Section 3.3.8."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [("(223|A)", (223, "A", None)), ("(225|C|dbSNP)", (225, "C", "dbSNP"))],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        (v,) = _entry(rf"\VariantSimple={item}").variant_simple
        assert (v.position, v.new_amino_acid, v.tag) == expected
        _assert_entry_round_trip(rf"\VariantSimple={item}")


class TestVariantComplexExamples:
    """Section 3.3.9 table."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [
            ("(100|100|)", (100, 100, "", None)),
            ("(100|100||10kexomes)", (100, 100, "", "10kexomes")),
            ("(100|102|)", (100, 102, "", None)),
            ("(100|100|APT)", (100, 100, "APT", None)),
            ("(100|102|KPA)", (100, 102, "KPA", None)),
            ("(100|101|P)", (100, 101, "P", None)),
        ],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        (v,) = _entry(rf"\VariantComplex={item}").variant_complex
        assert (v.start_pos, v.end_pos, v.new_sequence, v.tag) == expected
        _assert_entry_round_trip(rf"\VariantComplex={item}")


class TestModResUnimodExamples:
    """Section 3.3.10 table."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [
            ("(100|UNIMOD:21|Phospho)", ((100,), None)),
            ("(100,157,214|UNIMOD:21|Phospho)", ((100, 157, 214), None)),
            ("(100,157|UNIMOD:21|Phospho|invitro)", ((100, 157), "invitro")),
            ("(?|UNIMOD:21|Phospho)", (("?",), None)),
        ],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        (m,) = _entry(rf"\ModResUnimod={item}").mod_res_unimod
        assert (m.positions, m.accession, m.name, m.tag) == (expected[0], "UNIMOD:21", "Phospho", expected[1])
        _assert_entry_round_trip(rf"\ModResUnimod={item}")

    def test_repeated_unknown_position(self) -> None:
        # '"?" may be listed multiple times to denote multiple modifications'.
        (m,) = _entry(r"\ModResUnimod=(?,?|UNIMOD:21|Phospho)").mod_res_unimod
        assert m.positions == ("?", "?")


class TestModResPsiExamples:
    """Section 3.3.11 table."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [
            ("(100|MOD:00046|O-phospho-L-serine)", ((100,), None, None)),
            ("(12:100|MOD:00046|O-phospho-L-serine)", ((100,), None, 12)),
            ("(100,157|MOD:00046|O-phospho-L-serine)", ((100, 157), None, None)),
            ("(100,157,214|MOD:00046|O-phospho-L-serine|uncertain)", ((100, 157, 214), "uncertain", None)),
            ("(?|MOD:00046|O-phospho-L-serine)", (("?",), None, None)),
        ],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        (m,) = _entry(rf"\ModResPsi={item}").mod_res_psi
        assert (m.positions, m.accession, m.name, m.tag, m.annot_id) == (
            expected[0],
            "MOD:00046",
            "O-phospho-L-serine",
            expected[1],
            expected[2],
        )
        _assert_entry_round_trip(rf"\ModResPsi={item}")


class TestModResExamples:
    """Section 3.3.12 table."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [
            ("(100||N-linked (GlcNAc...))", ((100,), "", "N-linked (GlcNAc...)", None)),
            ("(100,178||N-linked (GlcNAc...)|invitro)", ((100, 178), "", "N-linked (GlcNAc...)", "invitro")),
            ("(100|CustomMod:22|Floxilation)", ((100,), "CustomMod:22", "Floxilation", None)),
            ("(100||Phosphorylation)", ((100,), "", "Phosphorylation", None)),
        ],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        (m,) = _entry(rf"\ModRes={item}").mod_res
        assert (m.positions, m.accession, m.name, m.tag) == expected
        _assert_entry_round_trip(rf"\ModRes={item}")


class TestProcessedExamples:
    """Section 3.3.13 table."""

    @pytest.mark.parametrize(
        ("item", "expected"),
        [
            ("(1|40|PEFF:0001021|signal peptide)", (1, 40, "PEFF:0001021", "signal peptide")),
            ("(41|390|PEFF:0001020|mature protein)", (41, 390, "PEFF:0001020", "mature protein")),
        ],
    )
    def test_legal(self, item: str, expected: tuple) -> None:
        # The spec's second example ends at 890 (TYRO3); 390 keeps it inside _SEQ.
        (p,) = _entry(rf"\Processed={item}").processed
        assert (p.start_pos, p.end_pos, p.accession, p.name, p.tag) == (*expected, None)
        _assert_entry_round_trip(rf"\Processed={item}")


# ---------------------------------------------------------------------------
# ILLEGAL examples: every one must be reported
# ---------------------------------------------------------------------------


class TestIllegalExamples:
    """Every value the spec tables mark ILLEGAL gives a PeffWarning or PeffParseError."""

    @pytest.mark.parametrize(
        ("desc", "match"),
        [
            # 3.3.9: "Not a legal VariantComplex. This MUST be encoded as a VariantSimple."
            (r"\VariantComplex=(100|100|A)", "VariantSimple"),
            # 3.3.9: "No regular expressions are allowed in this item."
            (r"\VariantComplex=(100|100|[AEQ]P)", "newSequence"),
            # 3.3.10: accession and name MUST be provided.
            (r"\ModResUnimod=(100||Phospho)", "accession must be provided"),
            (r"\ModResUnimod=(100|UNIMOD:21|)", "name must be provided"),
            # 3.3.11
            (r"\ModResPsi=(100||O-phospho-L-serine)", "accession must be provided"),
            (r"\ModResPsi=(100|MOD:00046|)", "name must be provided"),
            # 3.3.13: accession and name from the PEFF CV MUST be provided.
            (r"\Processed=(1|40||signal peptide)", "accession must be provided"),
            (r"\Processed=(1|40|PEFF:0001021|)", "name must be provided"),
        ],
    )
    def test_warns(self, desc: str, match: str) -> None:
        with pytest.warns(PeffWarning, match=match):
            read_peff(StringIO(f"{_HEADER}>sp:P1 {desc}\n{_SEQ}\n"))

    def test_mod_res_missing_accession_component(self) -> None:
        # 3.3.12: "(100|Phosphorylation) ILLEGAL ... skipping the second element is
        # not permitted."
        with pytest.raises(PeffParseError, match="ModRes"):
            read_peff(StringIO(f"{_HEADER}>sp:P1 \\ModRes=(100|Phosphorylation)\n{_SEQ}\n"))

    @pytest.mark.parametrize("new_aa", ["", " ", "-", "AC", "1"])
    def test_variant_simple_bad_amino_acid(self, new_aa: str) -> None:
        # 3.3.8: newAminoAcid "MUST be a valid amino acid code ... or an asterisk (*).
        # It MUST NOT be empty, or space, or any non-alphabetic character except asterisk."
        with pytest.warns(PeffWarning, match="newAminoAcid"):
            read_peff(StringIO(f"{_HEADER}>sp:P1 \\VariantSimple=(5|{new_aa})\n{_SEQ}\n"))

    @pytest.mark.parametrize("new_aa", ["A", "J", "X", "*", "u"])
    def test_variant_simple_good_amino_acid(self, new_aa: str) -> None:
        assert _entry(rf"\VariantSimple=(5|{new_aa})").variant_simple[0].new_amino_acid == new_aa


# ---------------------------------------------------------------------------
# 3.4.2 annotation identifiers: the insulin example
# ---------------------------------------------------------------------------


class TestInsulinExample:
    """Section 3.4.2: 11 insulin proteoforms via annotation identifiers."""

    def _read(self) -> tuple[FileHeader, list[SequenceEntry]]:
        with warnings.catch_warnings():
            warnings.simplefilter("error", PeffWarning)
            return read_peff(SPEC / "spec_3_4_2_insulin.peff")

    def test_header_flag(self) -> None:
        header, _ = self._read()
        assert header.databases[0].has_annotation_identifiers is True

    def test_annotation_ids_are_consecutive(self) -> None:
        (entry,) = self._read()[1]
        ids = [
            a.annot_id
            for group in (
                entry.mod_res_psi,
                entry.variant_simple,
                entry.processed,
                entry.disulfide_bond,
            )
            for a in group
        ]
        assert ids == list(range(84))

    def test_mod_res_psi(self) -> None:
        (entry,) = self._read()[1]
        assert entry.mod_res_psi[0].positions == (53,)
        assert entry.mod_res_psi[0].name == "N6-myristoyl-L-lysine"
        assert [(m.annot_id, m.positions[0]) for m in entry.mod_res_psi[1:]] == [
            (1, 31),
            (2, 96),
            (3, 43),
            (4, 109),
            (5, 95),
            (6, 100),
        ]

    def test_variant_simple(self) -> None:
        (entry,) = self._read()[1]
        assert len(entry.variant_simple) == 70
        assert (entry.variant_simple[0].position, entry.variant_simple[0].new_amino_acid) == (2, "T")
        assert (entry.variant_simple[-1].position, entry.variant_simple[-1].new_amino_acid) == (108, "C")

    def test_processed(self) -> None:
        (entry,) = self._read()[1]
        assert [(p.annot_id, p.start_pos, p.end_pos, p.accession, p.name) for p in entry.processed] == [
            (77, 1, 24, "PEFF:0001021", "signal peptide"),
            (78, 25, 54, "PEFF:0001020", "mature protein"),
            (79, 57, 87, "PEFF:0001034", "propeptide"),
            (80, 90, 110, "PEFF:0001020", "mature protein"),
        ]

    def test_disulfide_bonds_reference_mod_res_psi_ids(self) -> None:
        (entry,) = self._read()[1]
        assert [(d.annot_id, d.annot_id_refs, d.description) for d in entry.disulfide_bond] == [
            (81, (1, 2), "between chains"),
            (82, (3, 4), "between chains"),
            (83, (5, 6), "A chain only"),
        ]

    def test_proteoforms(self) -> None:
        (entry,) = self._read()[1]
        pfs = entry.proteoform
        assert [p.proteoform_id for p in pfs] == [f"NX_P01308-1-pf{i}" for i in range(1, 12)]
        assert [(r.start, r.end) for r in pfs[0].ranges] == [(1, 110)]
        assert pfs[0].annot_id_refs == ()
        assert pfs[0].name == "preproinsulin"
        assert pfs[2].annot_id_refs == (1, 2, 3, 4, 5, 6)
        assert pfs[7].annot_id_refs == (0, 1, 3)
        assert [(r.start, r.end) for r in pfs[10].ranges] == [(90, 110), (25, 54)]
        assert pfs[10].annot_id_refs == (81, 82, 83)
        assert pfs[10].name == "Insulin: chains A and B joined"

    def test_short_disulfide_example(self) -> None:
        # "\DisulfideBond=(15:1,2|between chains)"
        (bond,) = _entry(r"\DisulfideBond=(15:1,2|between chains)").disulfide_bond
        assert (bond.annot_id, bond.annot_id_refs, bond.description) == (15, (1, 2), "between chains")

    def test_writes_back(self) -> None:
        _assert_round_trip(*self._read())


# ---------------------------------------------------------------------------
# 3.3.5 optional tags and section 8 glossary forms
# ---------------------------------------------------------------------------


class TestOptionalTagExamples:
    """Section 3.3.5: example tags and the OptionalTagDef header key."""

    @pytest.mark.parametrize("tag", ["uncertain", "dbSNP", "in vitro", "[sample04][sample08]"])
    def test_tag_values(self, tag: str) -> None:
        assert _entry(rf"\VariantSimple=(5|A|{tag})").variant_simple[0].tag == tag
        _assert_entry_round_trip(rf"\VariantSimple=(5|A|{tag})")

    def test_no_trailing_pipe_without_tag(self) -> None:
        # "If such a tag is not provided, the trailing pipe character MUST NOT be written."
        header, entries = _read_strict(f"{_HEADER}>sp:P1 \\VariantSimple=(5|A)\n{_SEQ}\n")
        assert "\\VariantSimple=(5|A)\n" in _write(header, entries)

    def test_optional_tag_def(self) -> None:
        text = _HEADER.removesuffix("# //\n") + "# OptionalTagDef=dbSNP:variant from dbSNP\n# //\n"
        header = _read_header_strict(text)
        (otd,) = header.databases[0].optional_tag_defs
        assert (otd.tag, otd.description) == ("dbSNP", "variant from dbSNP")
        _assert_round_trip(header, [])


class TestGlossaryForms:
    """Section 8: a single item or component list may be written without parentheses."""

    def test_single_unparenthesized_item(self) -> None:
        # \Key=Component1|Component2
        (m,) = _entry(r"\ModResUnimod=100|UNIMOD:21|Phospho").mod_res_unimod
        assert (m.positions, m.accession, m.name) == ((100,), "UNIMOD:21", "Phospho")

    def test_parenthesized_items_with_tags(self) -> None:
        # \Key=(Component1|Component2|OptionalTag)(Component1|Component2|OptionalTag)
        entry = _entry(r"\VariantSimple=(5|A|t1)(6|C|t2)")
        assert [(v.position, v.new_amino_acid, v.tag) for v in entry.variant_simple] == [(5, "A", "t1"), (6, "C", "t2")]
