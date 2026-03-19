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
