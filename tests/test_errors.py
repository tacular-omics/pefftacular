"""Tests for error handling on malformed PEFF input."""

import warnings
from io import StringIO

import pytest

from pefftacular._parser import PeffReader
from pefftacular.errors import PeffError, PeffParseError, PeffWarning, PeffWriteError

# ---------------------------------------------------------------------------
# Minimal PEFF header shared across integer-parsing and length-mismatch tests
# ---------------------------------------------------------------------------

_HEADER = (
    "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
    "# NumberOfEntries=1\n# SequenceType=AA\n# //\n"
)


class TestMalformedHeader:
    def test_no_peff_line(self):
        with pytest.raises(PeffParseError, match="First line"):
            PeffReader(StringIO(">sp:X1\nACDEF\n")).header  # noqa: B018

    def test_garbage_first_line(self):
        with pytest.raises(PeffParseError, match="First line"):
            PeffReader(StringIO("garbage\n")).header  # noqa: B018

    def test_empty_file(self):
        with pytest.raises(PeffParseError, match="Empty file"):
            PeffReader(StringIO("")).header  # noqa: B018

    def test_whitespace_only(self):
        with pytest.raises(PeffParseError):
            PeffReader(StringIO("   \n  \n")).header  # noqa: B018


class TestMalformedEntry:
    def test_bad_description_line(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n>NOCOLON \\Length=3\nABC\n"
        )
        with pytest.raises(PeffParseError, match="Invalid description"):
            list(PeffReader(StringIO(data)))

    def test_variant_simple_too_few_fields(self):
        data = (
            "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n"
            "# NumberOfEntries=1\n# SequenceType=AA\n# //\n>t:X1 \\VariantSimple=(5)\nACDEF\n"
        )
        with pytest.raises(PeffParseError, match="VariantSimple"):
            list(PeffReader(StringIO(data)))


class TestIntegerParsingErrors:
    def test_number_of_entries_non_integer_raises(self) -> None:
        data = (
            "# PEFF 1.0\n"
            "# //\n"
            "# DbName=t\n"
            "# Prefix=t\n"
            "# DbVersion=1\n"
            "# DbSource=x\n"
            "# NumberOfEntries=not_a_number\n"
            "# SequenceType=AA\n"
            "# //\n"
        )
        with pytest.raises(PeffParseError, match="NumberOfEntries"):
            PeffReader(StringIO(data)).header  # noqa: B018

    def test_ncbi_tax_id_non_integer_raises(self) -> None:
        data = _HEADER + ">t:X1 \\NcbiTaxId=bad\nACDEF\n"
        with pytest.raises(PeffParseError, match="NcbiTaxId"):
            list(PeffReader(StringIO(data)))

    def test_length_non_integer_raises(self) -> None:
        data = _HEADER + ">t:X1 \\Length=abc\nACDEF\n"
        with pytest.raises(PeffParseError, match="Length"):
            list(PeffReader(StringIO(data)))

    def test_sv_non_integer_raises(self) -> None:
        data = _HEADER + ">t:X1 \\SV=x\nACDEF\n"
        with pytest.raises(PeffParseError, match="SV"):
            list(PeffReader(StringIO(data)))

    def test_ev_non_integer_raises(self) -> None:
        data = _HEADER + ">t:X1 \\EV=x\nACDEF\n"
        with pytest.raises(PeffParseError, match="EV"):
            list(PeffReader(StringIO(data)))

    def test_pe_non_integer_raises(self) -> None:
        data = _HEADER + ">t:X1 \\PE=x\nACDEF\n"
        with pytest.raises(PeffParseError, match="PE"):
            list(PeffReader(StringIO(data)))


class TestLengthMismatchWarning:
    def test_length_mismatch_emits_user_warning(self) -> None:
        # Length=5 but sequence "MKTLL" has 5 chars — use a longer sequence to force mismatch
        data = _HEADER + ">t:X1 \\Length=5\nMKTLLMKTLL\n"
        with pytest.warns(UserWarning, match="Length=5"):
            list(PeffReader(StringIO(data)))

    def test_length_match_no_warning(self) -> None:
        data = _HEADER + ">t:X1 \\Length=5\nMKTLL\n"
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            # Should not raise — sequence length matches Length annotation
            list(PeffReader(StringIO(data)))


class TestNumberOfEntriesMismatch:
    def test_mismatch_emits_warning(self) -> None:
        _header_n5 = _HEADER.replace("# NumberOfEntries=1\n", "# NumberOfEntries=5\n")
        data = _header_n5 + ">t:X1 \\Length=2\nAC\n>t:X2 \\Length=2\nAC\n"
        with pytest.warns(UserWarning, match="NumberOfEntries=5"):
            list(PeffReader(StringIO(data)))

    def test_match_emits_no_warning(self) -> None:
        data = _HEADER + ">t:X1 \\Length=2\nAC\n"
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            list(PeffReader(StringIO(data)))


class TestEntryLineNumberIsAbsolute:
    def test_entry_error_reports_absolute_line(self) -> None:
        # Header is 9 lines; bad SV is on the 10th absolute line.
        data = _HEADER + ">t:X1 \\SV=bad\nAC\n"
        with pytest.raises(PeffParseError) as exc:
            list(PeffReader(StringIO(data)))
        assert exc.value.line == 10


class TestAnnotationValidationWarnings:
    """Spec MUST-rules (sections 3.3.8-3.3.13) surface as warnings, not errors."""

    def test_position_out_of_range_warns(self) -> None:
        # Sequence "MKTLL" is 5 residues; position 99 is out of range.
        data = _HEADER + ">t:X1 \\VariantSimple=(99|A)\nMKTLL\n"
        with pytest.warns(UserWarning, match="out of range 1..5"):
            list(PeffReader(StringIO(data)))

    def test_empty_new_amino_acid_warns(self) -> None:
        data = _HEADER + ">t:X1 \\VariantSimple=(3|)\nMKTLL\n"
        with pytest.warns(UserWarning, match="newAminoAcid must not be empty"):
            list(PeffReader(StringIO(data)))

    def test_unimod_missing_accession_warns(self) -> None:
        data = _HEADER + ">t:X1 \\ModResUnimod=(2||Phospho)\nMKTLL\n"
        with pytest.warns(UserWarning, match="ModResUnimod accession must be provided"):
            list(PeffReader(StringIO(data)))

    def test_valid_annotations_emit_no_warning(self) -> None:
        data = _HEADER + ">t:X1 \\VariantSimple=(3|A) \\ModResUnimod=(2|UNIMOD:21|Phospho)\nMKTLL\n"
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            list(PeffReader(StringIO(data)))

    def test_disulfide_refs_not_range_checked(self) -> None:
        # DisulfideBond values are annotation-ID references, not residue
        # positions, so high values must NOT trigger a range warning.
        data = _HEADER + ">t:X1 \\DisulfideBond=(100,200|refs)\nMKTLL\n"
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            list(PeffReader(StringIO(data)))


class TestPeffParseErrorAttributes:
    def test_line_number_in_message(self):
        err = PeffParseError("bad thing", line=42)
        assert "Line 42" in str(err)
        assert err.line == 42

    def test_no_line_number(self):
        err = PeffParseError("bad thing")
        assert "Line" not in str(err)
        assert err.line is None

    def test_context_stored(self):
        err = PeffParseError("bad", context="the offending line")
        assert err.context == "the offending line"


class TestErrorHierarchyAndHints:
    """PeffError base class, hints, and exception notes for coding agents."""

    def test_parse_error_is_peff_error_and_value_error(self):
        err = PeffParseError("boom")
        assert isinstance(err, PeffError)
        assert isinstance(err, ValueError)

    def test_write_error_is_peff_error(self):
        err = PeffWriteError("boom")
        assert isinstance(err, PeffError)
        assert isinstance(err, ValueError)

    def test_hint_stored_and_attached_as_note(self):
        err = PeffParseError("bad", context="X", hint="do Y instead")
        assert err.hint == "do Y instead"
        # Notes surface in tracebacks (PEP 678) but not in str().
        assert any("do Y instead" in note for note in err.__notes__)
        assert any("offending text" in note for note in err.__notes__)
        assert "do Y instead" not in str(err)

    def test_parse_failure_carries_actionable_hint(self):
        data = _HEADER + ">t:X1 \\VariantSimple=(5)\nMKTLL\n"
        with pytest.raises(PeffParseError) as exc:
            list(PeffReader(StringIO(data)))
        assert exc.value.hint is not None
        assert "position|newAminoAcid" in exc.value.hint

    def test_write_error_carries_hint(self):
        with pytest.raises(PeffWriteError) as exc:
            from pefftacular import FileHeader, SequenceEntry, write_peff

            write_peff(FileHeader(peff_version="1.0"), [SequenceEntry("p", "id", "")], StringIO())
        assert exc.value.hint is not None


class TestPeffWarningCategory:
    """Spec-violation warnings use the PeffWarning category (a UserWarning subclass)."""

    def test_out_of_range_uses_peff_warning(self):
        data = _HEADER + ">t:X1 \\VariantSimple=(99|A)\nMKTLL\n"
        with pytest.warns(PeffWarning, match="out of range"):
            list(PeffReader(StringIO(data)))

    def test_peff_warning_is_user_warning(self):
        assert issubclass(PeffWarning, UserWarning)

    def test_can_escalate_only_peff_warnings(self):
        data = _HEADER + ">t:X1 \\VariantSimple=(99|A)\nMKTLL\n"
        with warnings.catch_warnings():  # noqa: SIM117
            warnings.simplefilter("error", PeffWarning)
            with pytest.raises(PeffWarning):
                list(PeffReader(StringIO(data)))
