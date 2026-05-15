"""Tests for PEFF writer."""

import io
from io import StringIO

import pytest

from pefftacular._models import (
    CustomKeyDef,
    DatabaseHeader,
    DisulfideBond,
    FileHeader,
    ModResUnimod,
    OptionalTagDef,
    Processed,
    Proteoform,
    SequenceEntry,
    SequenceRange,
    VariantComplex,
    VariantSimple,
)
from pefftacular._writer import write_peff
from pefftacular.errors import PeffWriteError


def _make_minimal_header() -> FileHeader:
    return FileHeader(
        peff_version="1.0",
        databases=(
            DatabaseHeader(
                prefix="sp",
                db_name="testdb",
                db_version="2024-01",
                db_sources=("http://example.com",),
                number_of_entries=1,
                sequence_type="AA",
            ),
        ),
    )


class TestHeaderWriting:
    def test_version_line(self):
        buf = StringIO()
        write_peff(_make_minimal_header(), [], buf)
        lines = buf.getvalue().splitlines()
        assert lines[0] == "# PEFF 1.0"

    def test_general_comments(self):
        header = FileHeader(
            peff_version="1.0",
            general_comments=("Comment one", "Comment two"),
            databases=(
                DatabaseHeader(prefix="sp", db_version="1", db_sources=("x",), number_of_entries=0, sequence_type="AA"),
            ),
        )
        buf = StringIO()
        write_peff(header, [], buf)
        text = buf.getvalue()
        assert "# GeneralComment=Comment one\n" in text
        assert "# GeneralComment=Comment two\n" in text

    def test_database_block(self):
        buf = StringIO()
        write_peff(_make_minimal_header(), [], buf)
        text = buf.getvalue()
        assert "# DbName=testdb\n" in text
        assert "# Prefix=sp\n" in text
        assert "# DbVersion=2024-01\n" in text
        assert "# NumberOfEntries=1\n" in text
        assert "# SequenceType=AA\n" in text

    def test_custom_key_def(self):
        header = FileHeader(
            peff_version="1.0",
            databases=(
                DatabaseHeader(
                    prefix="sp",
                    db_version="1",
                    db_sources=("x",),
                    number_of_entries=0,
                    sequence_type="AA",
                    custom_key_defs=(
                        CustomKeyDef(
                            key_name="MyKey",
                            description="A custom key",
                            field_names=("a", "b"),
                            field_types=("string", "integer"),
                        ),
                    ),
                ),
            ),
        )
        buf = StringIO()
        write_peff(header, [], buf)
        text = buf.getvalue()
        assert "CustomKeyDef=" in text
        assert "KeyName=MyKey" in text


class TestEntryWriting:
    def test_simple_entry(self):
        entry = SequenceEntry(prefix="sp", db_unique_id="P12345", sequence="ACDEF", pname="Albumin", length=5)
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert ">sp:P12345 " in text
        assert "\\PName=Albumin" in text
        assert "\\Length=5" in text
        assert "ACDEF" in text

    def test_sequence_wrapping(self):
        seq = "A" * 130
        entry = SequenceEntry(prefix="sp", db_unique_id="X1", sequence=seq, length=130)
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        lines = buf.getvalue().splitlines()
        # Find sequence lines (after description)
        seq_lines = [ln for ln in lines if ln and not ln.startswith("#") and not ln.startswith(">")]
        assert len(seq_lines) == 3  # 60 + 60 + 10
        assert len(seq_lines[0]) == 60
        assert len(seq_lines[1]) == 60
        assert len(seq_lines[2]) == 10

    def test_variant_simple_serialization(self):
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="X1",
            sequence="ACDEF",
            variant_simple=(
                VariantSimple(position=5, new_amino_acid="A"),
                VariantSimple(position=10, new_amino_acid="G", tag="rs123"),
            ),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\VariantSimple=(5|A)(10|G|rs123)" in text

    def test_variant_complex_serialization(self):
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="X1",
            sequence="ACDEF",
            variant_complex=(VariantComplex(start_pos=3, end_pos=7, new_sequence="AAA", tag="t1"),),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\VariantComplex=(3|7|AAA|t1)" in text

    def test_mod_res_unimod_positions(self):
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="X1",
            sequence="ACDEF",
            mod_res_unimod=(ModResUnimod(positions=(42, 57), accession="UNIMOD:21", name="Phospho"),),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\ModResUnimod=(42,57|UNIMOD:21|Phospho)" in text

    def test_processed_serialization(self):
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="X1",
            sequence="ACDEF",
            processed=(Processed(start_pos=1, end_pos=20, accession="PEFF:0001", name="mature protein"),),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\Processed=(1|20|PEFF:0001|mature protein)" in text

    def test_decoy_bool(self):
        entry = SequenceEntry(prefix="sp", db_unique_id="X1", sequence="AC", decoy=True)
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        assert "\\Decoy=true" in buf.getvalue()

    def test_none_fields_skipped(self):
        entry = SequenceEntry(prefix="sp", db_unique_id="X1", sequence="AC")
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\PName" not in text
        assert "\\GName" not in text
        assert "\\Length" not in text

    def test_canonical_key_order(self):
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="X1",
            sequence="ACDEF",
            pname="Test",
            length=5,
            gname="TST",
            ncbi_tax_id=9606,
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        desc_line = [ln for ln in buf.getvalue().splitlines() if ln.startswith(">")][0]
        # Length should come before PName
        assert desc_line.index("\\Length") < desc_line.index("\\PName")
        # PName before GName
        assert desc_line.index("\\PName") < desc_line.index("\\GName")
        # GName before NcbiTaxId
        assert desc_line.index("\\GName") < desc_line.index("\\NcbiTaxId")


class TestHeaderOptionalFields:
    def test_all_optional_database_fields_serialized(self) -> None:
        header = FileHeader(
            peff_version="1.0",
            databases=(
                DatabaseHeader(
                    prefix="sp",
                    db_name="testdb",
                    db_description="A description",
                    db_version="2024-01",
                    db_date="20240101",
                    db_sources=("x",),
                    number_of_entries=0,
                    sequence_type="AA",
                    decoy=False,
                    conversion="some-conversion",
                    has_annotation_identifiers=True,
                    proteoform_db=True,
                    optional_tag_defs=(OptionalTagDef(tag="t1", description="desc1"),),
                    extra={"CustomDbKey": "custom_value"},
                ),
            ),
        )
        buf = StringIO()
        write_peff(header, [], buf)
        text = buf.getvalue()
        assert "# DbDescription=A description\n" in text
        assert "# DbDate=20240101\n" in text
        assert "# Decoy=false\n" in text
        assert "# Conversion=some-conversion\n" in text
        assert "# HasAnnotationIdentifiers=true\n" in text
        assert "# ProteoformDb=true\n" in text
        assert "# OptionalTagDef=t1:desc1\n" in text
        assert "# CustomDbKey=custom_value\n" in text


class TestEntryOptionalFields:
    def test_comment_serialized(self) -> None:
        entry = SequenceEntry(
            prefix="sp", db_unique_id="P1", sequence="AC", comment="An entry comment"
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        assert "\\Comment=An entry comment" in buf.getvalue()

    def test_disulfide_bond_serialized(self) -> None:
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="P1",
            sequence="ACDEF",
            disulfide_bond=(
                DisulfideBond(positions=(1, 2), description="between chains", annot_id=81),
                DisulfideBond(positions=(3, 4)),
            ),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\DisulfideBond=(81:1,2|between chains)(3,4)" in text

    def test_proteoform_serialized(self) -> None:
        entry = SequenceEntry(
            prefix="sp",
            db_unique_id="P1",
            sequence="ACDEF",
            proteoform=(
                Proteoform(
                    proteoform_id="pf1",
                    ranges=(SequenceRange(start=1, end=10),),
                    annot_id_refs=(1, 2),
                    name="some form",
                    annot_id=100,
                ),
                Proteoform(
                    proteoform_id="pf2",
                    ranges=(SequenceRange(start=5, end=20),),
                    annot_id_refs=(),
                ),
            ),
        )
        buf = StringIO()
        write_peff(_make_minimal_header(), [entry], buf)
        text = buf.getvalue()
        assert "\\Proteoform=(100:pf1|1-10|1,2|some form)(pf2|5-20|)" in text


class TestWritePeffValidation:
    def _valid_entry(self) -> SequenceEntry:
        return SequenceEntry(prefix="sp", db_unique_id="P00001", sequence="ACDEF")

    def test_none_header_raises(self) -> None:
        with pytest.raises(PeffWriteError, match="header"):
            write_peff(None, [], io.StringIO())  # type: ignore[arg-type]

    def test_empty_prefix_raises(self) -> None:
        entry = SequenceEntry(prefix="", db_unique_id="P00001", sequence="ACDEF")
        with pytest.raises(PeffWriteError, match="prefix"):
            write_peff(_make_minimal_header(), [entry], io.StringIO())

    def test_empty_db_unique_id_raises(self) -> None:
        entry = SequenceEntry(prefix="sp", db_unique_id="", sequence="ACDEF")
        with pytest.raises(PeffWriteError, match="db_unique_id"):
            write_peff(_make_minimal_header(), [entry], io.StringIO())

    def test_empty_sequence_raises(self) -> None:
        entry = SequenceEntry(prefix="sp", db_unique_id="P00001", sequence="")
        with pytest.raises(PeffWriteError, match="sequence"):
            write_peff(_make_minimal_header(), [entry], io.StringIO())
