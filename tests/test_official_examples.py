"""Tests for the official PEFF example fixtures."""

import warnings
from pathlib import Path

import pytest

from pefftacular._parser import PeffReader, read_peff
from pefftacular.errors import PeffParseError

FIXTURES = Path(__file__).parent / "fixtures"


class TestMinimalValid:
    """PEFF_Minimal_Valid.peff — single-entry minimal file."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "PEFF_Minimal_Valid.peff") as reader:
            list(reader)

    def test_entry_count(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert len(entries) == 1

    def test_entry_prefix(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert entries[0].prefix == "sp"

    def test_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert entries[0].db_unique_id == "Q9Y2X3"

    def test_entry_length(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert entries[0].length == 1

    def test_entry_sequence(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert entries[0].sequence == "M"

    def test_peff_version(self) -> None:
        header, _ = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert header.peff_version == "1.0"

    def test_database_name(self) -> None:
        header, _ = read_peff(FIXTURES / "PEFF_Minimal_Valid.peff")
        assert header.databases[0].db_name == "Minimal Test example PEFF_Minimal_Valid.peff"


class TestTinyValid:
    """PEFF_Tiny_Valid.peff — three entries across three databases."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "PEFF_Tiny_Valid.peff") as reader:
            list(reader)

    def test_entry_count(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert len(entries) == 3

    def test_database_count(self) -> None:
        header, _ = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert len(header.databases) == 3

    def test_database_prefixes(self) -> None:
        header, _ = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        prefixes = [db.prefix for db in header.databases]
        assert prefixes == ["sp", "nxp", "nr"]

    def test_entry_prefixes_match_databases(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[0].prefix == "sp"
        assert entries[1].prefix == "nxp"
        assert entries[2].prefix == "nr"

    def test_sp_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[0].db_unique_id == "Q9Y2X3"

    def test_nxp_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[1].db_unique_id == "NX_P11171-2"

    def test_nxp_entry_gene_name(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[1].gname == "EPB41"

    def test_nxp_entry_has_mod_res_psi(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert len(entries[1].mod_res_psi) > 0

    def test_nr_entry_has_mod_res_unimod(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert len(entries[2].mod_res_unimod) > 0

    def test_sp_entry_ncbi_tax_id(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[0].ncbi_tax_id == 9606

    def test_sp_entry_length(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
        assert entries[0].length == 528


class TestSmallTestDbPeff10:
    """SmallTestDB-PEFF1.0.peff — 29 entries across 2 databases."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "SmallTestDB-PEFF1.0.peff") as reader:
            list(reader)

    def test_entry_count_matches_header_declarations(self) -> None:
        header, entries = read_peff(FIXTURES / "SmallTestDB-PEFF1.0.peff")
        declared = sum(db.number_of_entries for db in header.databases if db.number_of_entries is not None)
        assert len(entries) == declared

    def test_two_databases(self) -> None:
        header, _ = read_peff(FIXTURES / "SmallTestDB-PEFF1.0.peff")
        assert len(header.databases) == 2

    def test_sp_database_declares_28_entries(self) -> None:
        header, _ = read_peff(FIXTURES / "SmallTestDB-PEFF1.0.peff")
        sp_db = next(db for db in header.databases if db.prefix == "sp")
        assert sp_db.number_of_entries == 28

    def test_peff_version(self) -> None:
        header, _ = read_peff(FIXTURES / "SmallTestDB-PEFF1.0.peff")
        assert header.peff_version == "1.0"


class TestUniProtExport3Prot:
    """UniProtExport_3prot.peff — 3 entries from a UniProt export."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "UniProtExport_3prot.peff") as reader:
            list(reader)

    def test_entry_count(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert len(entries) == 3

    def test_db_name_is_uniprot(self) -> None:
        header, _ = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert header.databases[0].db_name == "UniProt"

    def test_first_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert entries[0].db_unique_id == "P21802-23"

    def test_first_entry_prefix(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert entries[0].prefix == "up"

    def test_first_entry_length(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert entries[0].length == 709

    def test_first_entry_has_variant_simple(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert len(entries[0].variant_simple) > 0

    def test_third_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "UniProtExport_3prot.peff")
        assert entries[2].db_unique_id == "A0A1B0GTB3"


class TestAnnotIdInsulinValid:
    """PEFF_AnnotID_Insulin_Valid.peff — insulin with annotation identifiers."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff") as reader:
            list(reader)

    def test_entry_count(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries) == 1

    def test_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert entries[0].db_unique_id == "NX_P01308-1"

    def test_has_mod_res_psi(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries[0].mod_res_psi) == 7

    def test_has_variant_simple(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries[0].variant_simple) == 70

    def test_has_processed(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries[0].processed) == 4

    def test_disulfide_bond_parsed(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries[0].disulfide_bond) == 3

    def test_proteoform_parsed(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert len(entries[0].proteoform) == 11

    def test_gene_name(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        assert entries[0].gname == "INS"

    def test_sequence_length(self) -> None:
        _, entries = read_peff(FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff")
        e = entries[0]
        assert len(e.sequence) == e.length


class TestProteoformEnst:
    """proteoform_ENST00000000412.peff — proteogenomics file with ProteoformDb header."""

    def test_parses_without_error(self) -> None:
        with PeffReader(FIXTURES / "proteoform_ENST00000000412.peff") as reader:
            list(reader)

    def test_entry_count(self) -> None:
        _, entries = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        assert len(entries) == 4

    def test_proteoform_db_header_flag(self) -> None:
        header, _ = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        assert header.databases[0].proteoform_db is True

    def test_first_entry_id(self) -> None:
        _, entries = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        assert entries[0].db_unique_id == "ENST00000000412-1"

    def test_entries_have_proteoform_annotation(self) -> None:
        _, entries = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        for entry in entries:
            assert len(entry.proteoform) > 0

    def test_first_entry_has_variant_simple(self) -> None:
        _, entries = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        assert len(entries[0].variant_simple) > 0

    def test_first_entry_has_variant_complex(self) -> None:
        _, entries = read_peff(FIXTURES / "proteoform_ENST00000000412.peff")
        assert len(entries[0].variant_complex) > 0


class TestMinimalInvalid1:
    """PEFF_Minimal_INValid1.peff — version 0.0, missing mandatory keys, bogus tag."""

    def test_emits_version_warning(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            read_peff(FIXTURES / "PEFF_Minimal_INValid1.peff")
        messages = [str(w.message) for w in caught]
        assert any("0.0" in m for m in messages)

    def test_emits_missing_key_warnings(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            read_peff(FIXTURES / "PEFF_Minimal_INValid1.peff")
        messages = [str(w.message) for w in caught]
        assert any("missing mandatory key" in m for m in messages)

    def test_still_produces_two_entries(self) -> None:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, entries = read_peff(FIXTURES / "PEFF_Minimal_INValid1.peff")
        assert len(entries) == 2

    def test_bogus_color_tag_lands_in_extra(self) -> None:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, entries = read_peff(FIXTURES / "PEFF_Minimal_INValid1.peff")
        assert entries[0].extra.get("Color") == "Blue"


class TestTinyInvalid1:
    """PEFF_Tiny_INValid1.peff — version 1.99, BogusTag on first entry, malformed Processed."""

    def test_emits_version_warning(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(PeffParseError):
                read_peff(FIXTURES / "PEFF_Tiny_INValid1.peff")
        messages = [str(w.message) for w in caught]
        assert any("1.99" in m for m in messages)

    def test_raises_on_malformed_processed(self) -> None:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            with pytest.raises(PeffParseError, match="Processed"):
                read_peff(FIXTURES / "PEFF_Tiny_INValid1.peff")

    def test_header_parses_despite_invalid_version(self) -> None:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            with PeffReader(FIXTURES / "PEFF_Tiny_INValid1.peff") as reader:
                header = reader.header
        assert header.peff_version == "1.99"


class TestSmallTestDbPeff09:
    """SmallTestDB-PEFF0.9.peff — legacy format with no version number on header line."""

    def test_raises_parse_error_for_missing_version(self) -> None:
        with pytest.raises(PeffParseError, match="First line"):
            read_peff(FIXTURES / "SmallTestDB-PEFF0.9.peff")
