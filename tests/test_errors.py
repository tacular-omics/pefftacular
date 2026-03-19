"""Tests for error handling on malformed PEFF input."""

import warnings
from io import StringIO

import pytest

from pefftacular._parser import PeffReader
from pefftacular.errors import PeffParseError

# ---------------------------------------------------------------------------
# Minimal PEFF header shared across integer-parsing and length-mismatch tests
# ---------------------------------------------------------------------------

_HEADER = (
    "# PEFF 1.0\n"
    "# //\n"
    "# Prefix=t\n"
    "# DbVersion=1\n"
    "# DbSource=x\n"
    "# NumberOfEntries=1\n"
    "# SequenceType=AA\n"
    "# //\n"
)


class TestMalformedHeader:
    def test_no_peff_line(self):
        with pytest.raises(PeffParseError, match="First line"):
            PeffReader(StringIO(">sp:X1\nACDEF\n")).header

    def test_garbage_first_line(self):
        with pytest.raises(PeffParseError, match="First line"):
            PeffReader(StringIO("garbage\n")).header

    def test_empty_file(self):
        with pytest.raises(PeffParseError, match="Empty file"):
            PeffReader(StringIO("")).header

    def test_whitespace_only(self):
        with pytest.raises(PeffParseError):
            PeffReader(StringIO("   \n  \n")).header


class TestMalformedEntry:
    def test_bad_description_line(self):
        data = "# PEFF 1.0\n# //\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n# //\n>NOCOLON \\Length=3\nABC\n"
        with pytest.raises(PeffParseError, match="Invalid description"):
            list(PeffReader(StringIO(data)))

    def test_variant_simple_too_few_fields(self):
        data = "# PEFF 1.0\n# //\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n# //\n>t:X1 \\VariantSimple=(5)\nACDEF\n"
        with pytest.raises(PeffParseError, match="VariantSimple"):
            list(PeffReader(StringIO(data)))


class TestIntegerParsingErrors:
    def test_number_of_entries_non_integer_raises(self) -> None:
        data = (
            "# PEFF 1.0\n"
            "# //\n"
            "# Prefix=t\n"
            "# DbVersion=1\n"
            "# DbSource=x\n"
            "# NumberOfEntries=not_a_number\n"
            "# SequenceType=AA\n"
            "# //\n"
        )
        with pytest.raises(PeffParseError, match="NumberOfEntries"):
            PeffReader(StringIO(data)).header

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
