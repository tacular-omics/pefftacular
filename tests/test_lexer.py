"""Tests for the depth-tracking lexer."""

import pytest

from pefftacular._lexer import split_description_keys, split_fields, split_items
from pefftacular.errors import PeffParseError

# ---------------------------------------------------------------------------
# split_items
# ---------------------------------------------------------------------------


class TestSplitItems:
    def test_single_unparenthesized(self):
        assert split_items("110") == ["110"]

    def test_empty_string(self):
        assert split_items("") == [""]

    def test_single_parenthesized(self):
        assert split_items("(42|A)") == ["42|A"]

    def test_multiple_items(self):
        assert split_items("(42|A)(57|G)") == ["42|A", "57|G"]

    def test_nested_parens(self):
        result = split_items("(380||N-linked (GlcNAc...))(42|K)")
        assert result == ["380||N-linked (GlcNAc...)", "42|K"]

    def test_deeply_nested(self):
        result = split_items("(1|foo (bar (baz)))(2|x)")
        assert result == ["1|foo (bar (baz))", "2|x"]

    def test_mismatched_open_paren(self):
        with pytest.raises(PeffParseError, match="Unclosed"):
            split_items("(42|A")

    def test_mismatched_close_paren(self):
        with pytest.raises(PeffParseError, match="Unexpected"):
            split_items(")")


# ---------------------------------------------------------------------------
# split_fields
# ---------------------------------------------------------------------------


class TestSplitFields:
    def test_simple(self):
        assert split_fields("42|A") == ["42", "A"]

    def test_empty_middle(self):
        assert split_fields("380||N-linked (GlcNAc...)") == ["380", "", "N-linked (GlcNAc...)"]

    def test_no_pipe(self):
        assert split_fields("hello") == ["hello"]

    def test_pipe_inside_parens_ignored(self):
        assert split_fields("a|(b|c)|d") == ["a", "(b|c)", "d"]

    def test_four_fields(self):
        assert split_fields("42|UNIMOD:21|Phospho|mytag") == ["42", "UNIMOD:21", "Phospho", "mytag"]

    def test_trailing_pipe(self):
        assert split_fields("a|b|") == ["a", "b", ""]


# ---------------------------------------------------------------------------
# split_description_keys
# ---------------------------------------------------------------------------


class TestSplitDescriptionKeys:
    def test_simple_keys(self):
        result = split_description_keys(r"\PName=Albumin \GName=ALB \Length=10")
        assert result == {"PName": "Albumin", "GName": "ALB", "Length": "10"}

    def test_value_with_spaces(self):
        result = split_description_keys(r"\PName=Tyrosine-protein kinase receptor TYRO3 \GName=TYRO3")
        assert result == {"PName": "Tyrosine-protein kinase receptor TYRO3", "GName": "TYRO3"}

    def test_value_with_parens(self):
        result = split_description_keys(
            r"\ModResPsi=(681|MOD:00048|O4'-phospho-L-tyrosine)(686|MOD:00048|O4'-phospho-L-tyrosine) \Length=10"
        )
        assert result["ModResPsi"] == "(681|MOD:00048|O4'-phospho-L-tyrosine)(686|MOD:00048|O4'-phospho-L-tyrosine)"
        assert result["Length"] == "10"

    def test_empty_string(self):
        assert split_description_keys("") == {}

    def test_key_without_value(self):
        result = split_description_keys(r"\Decoy")
        assert result == {"Decoy": ""}

    def test_nested_parens_in_value(self):
        result = split_description_keys(r"\ModRes=(380||N-linked (GlcNAc...)) \Length=5")
        assert result["ModRes"] == "(380||N-linked (GlcNAc...))"
        assert result["Length"] == "5"

    def test_escaped_parens_do_not_change_depth(self):
        # An escaped unbalanced paren must not swallow the following key.
        result = split_description_keys(r"\ModRes=(1|X|odd \( name) \Length=5")
        assert result == {"ModRes": r"(1|X|odd \( name)", "Length": "5"}

    def test_escaped_close_paren_does_not_change_depth(self):
        result = split_description_keys(r"\ModRes=(1|X|odd \) name) \Length=5")
        assert result == {"ModRes": r"(1|X|odd \) name)", "Length": "5"}

    def test_escaped_backslash_after_space_is_not_a_key(self):
        result = split_description_keys(r"\PName=a \\ID=b \Length=5")
        assert result == {"PName": r"a \\ID=b", "Length": "5"}


class TestRepeatedKeysInLexer:
    """``split_description_keys`` keeps last-wins and reports what it dropped."""

    def test_mapping_keeps_the_last_value(self):
        assert split_description_keys(r"\PName=a \PName=b") == {"PName": "b"}

    def test_repeats_collects_every_earlier_value_in_file_order(self):
        repeats: list[tuple[str, str]] = []
        keys = split_description_keys(r"\ModRes=(1||a) \PName=x \ModRes=(2||b) \ModRes=(3||c)", repeats)
        assert keys == {"ModRes": "(3||c)", "PName": "x"}
        assert repeats == [("ModRes", "(1||a)"), ("ModRes", "(2||b)")]

    def test_repeats_is_empty_without_a_repeat(self):
        repeats: list[tuple[str, str]] = []
        split_description_keys(r"\PName=a \Length=5", repeats)
        assert repeats == []

    def test_escaped_path_reports_repeats_too(self):
        # An escape in the value forces the scanning path, not the fast path.
        repeats: list[tuple[str, str]] = []
        keys = split_description_keys(r"\ModRes=(1||a\|b) \ModRes=(2||c)", repeats)
        assert keys == {"ModRes": "(2||c)"}
        assert repeats == [("ModRes", r"(1||a\|b)")]

    def test_valueless_key_repeated(self):
        repeats: list[tuple[str, str]] = []
        assert split_description_keys(r"\Decoy \Decoy", repeats) == {"Decoy": ""}
        assert repeats == [("Decoy", "")]
