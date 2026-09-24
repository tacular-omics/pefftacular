"""The 1.1 lexer fast paths give the same result as the reference implementations.

The ``_ref_*`` functions below are the 1.0 lexer, copied verbatim (docstrings
dropped). Each fast path must return the same value, or raise the same error, for
any input, in particular inputs dense in ``\\ ( ) | " =`` and spaces.
"""

from __future__ import annotations

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pefftacular._lexer import _split_fields_escaped, split_description_keys, split_items
from pefftacular.errors import PeffParseError

# ---------------------------------------------------------------------------
# 1.0 reference implementations
# ---------------------------------------------------------------------------

_REF_ESCAPABLE = ("|", "(", ")", "\\")


def _ref_unescape_component(s: str) -> str:
    if "\\" not in s:
        return s
    out: list[str] = []
    i = 0
    length = len(s)
    while i < length:
        ch = s[i]
        if ch == "\\" and i + 1 < length and s[i + 1] in _REF_ESCAPABLE:
            out.append(s[i + 1])
            i += 2
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def _ref_has_unescaped(s: str, target: str) -> bool:
    i = 0
    length = len(s)
    while i < length:
        if s[i] == "\\" and i + 1 < length:
            i += 2
            continue
        if s[i] == target:
            return True
        i += 1
    return False


def _ref_split_items(raw: str, *, quotes: bool = False) -> list[str]:
    if not raw:
        return [raw]
    if not raw.startswith("("):
        if _ref_has_unescaped(raw, ")"):
            raise PeffParseError(
                "Unexpected ')' in value",
                context=raw,
                hint=r"Escape a literal ')' as '\)', or wrap multi-item values in matching parentheses",
            )
        return [raw]

    items: list[str] = []
    depth = 0
    item_start = -1
    in_quote = False
    i = 0
    length = len(raw)

    while i < length:
        ch = raw[i]
        if ch == "\\" and i + 1 < length:
            i += 2
            continue
        if quotes and ch == '"':
            in_quote = not in_quote
        elif in_quote:
            pass
        elif ch == "(":
            if depth == 0:
                item_start = i + 1
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                items.append(raw[item_start:i])
            elif depth < 0:
                raise PeffParseError(
                    f"Unexpected ')' at position {i}",
                    context=raw,
                    hint=r"Every ')' needs a matching '('; escape a literal one as '\)'",
                )
        i += 1

    if depth != 0:
        raise PeffParseError(
            "Unclosed '(' in value",
            context=raw,
            hint=r"Every '(' needs a matching ')'; escape a literal one as '\('",
        )

    return items


def _ref_split_fields_escaped(item: str) -> list[str]:
    fields: list[str] = []
    depth = 0
    start = 0
    i = 0
    length = len(item)

    while i < length:
        ch = item[i]
        if ch == "\\" and i + 1 < length:
            i += 2
            continue
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "|" and depth == 0:
            fields.append(_ref_unescape_component(item[start:i]))
            start = i + 1
        i += 1

    fields.append(_ref_unescape_component(item[start:]))
    return fields


def _ref_split_description_keys(rest: str) -> dict[str, str]:
    if not rest:
        return {}

    keys: dict[str, str] = {}
    depth = 0
    current_start: int | None = None  # index of the '\' that starts the current token
    i = 0
    length = len(rest)

    while i < length:
        ch = rest[i]
        if ch == "\\" and i + 1 < length and rest[i + 1] in _REF_ESCAPABLE:
            i += 2
            continue
        match ch:
            case "(":
                depth += 1
            case ")":
                depth -= 1
            case "\\" if depth == 0 and (i == 0 or rest[i - 1] == " "):
                # Close previous key-value if one is open
                if current_start is not None:
                    _ref_store_key_value(keys, rest[current_start:i].rstrip())
                current_start = i
        i += 1

    # Store last token
    if current_start is not None:
        _ref_store_key_value(keys, rest[current_start:])

    return keys


def _ref_store_key_value(keys: dict[str, str], token: str) -> None:
    # Strip leading backslash
    token = token.lstrip("\\")
    eq_idx = token.find("=")
    if eq_idx == -1:
        # Key with no value (e.g. ``\Decoy``)
        keys[token] = ""
    else:
        keys[token[:eq_idx]] = token[eq_idx + 1 :]


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------


def _outcome(fn, *args, **kwargs):  # type: ignore[no-untyped-def]
    try:
        return fn(*args, **kwargs)
    except PeffParseError as err:
        return ("error", str(err), err.context, err.hint)


_alphabet = st.sampled_from(list('\\()|" =ab1:?') + ["  ", "\\K", "(A|B)", "\\("])
_lexer_text = st.lists(_alphabet, max_size=30).map("".join) | st.text(max_size=40)


@given(_lexer_text, st.booleans())
@settings(max_examples=3000)
def test_split_items_matches_reference(raw: str, quotes: bool) -> None:
    assert _outcome(split_items, raw, quotes=quotes) == _outcome(_ref_split_items, raw, quotes=quotes)


@given(_lexer_text)
@settings(max_examples=3000)
def test_split_fields_escaped_matches_reference(item: str) -> None:
    assert _split_fields_escaped(item) == _ref_split_fields_escaped(item)


@given(_lexer_text)
@settings(max_examples=3000)
def test_split_description_keys_matches_reference(rest: str) -> None:
    assert split_description_keys(rest) == _ref_split_description_keys(rest)


@pytest.mark.parametrize(
    "rest",
    [
        "",
        "\\ID=P1 \\PName=Insulin \\Length=110",
        "\\VariantSimple=(12|A)(13|C|tag) \\ModResPsi=(5|MOD:00046|O-phospho-L-serine)",
        "\\PName=A (B \\C) \\GName=G",
        "\\PName=a\\(b \\GName=c",
        "\\Decoy \\ID=x",
        "\\PName=x  \\GName=y ",
        "text before \\ID=x",
        "\\PName=)) \\ID=x",
    ],
)
def test_split_description_keys_examples(rest: str) -> None:
    assert split_description_keys(rest) == _ref_split_description_keys(rest)


@pytest.mark.parametrize(
    "raw",
    ["", "110", "(A|B)(C|D)", "()", "(A)(B", "A)", "a\\)", "(A(B)C)", '("a(")', "(a\\)b)", "(x) (y)"],
)
@pytest.mark.parametrize("quotes", [False, True])
def test_split_items_examples(raw: str, quotes: bool) -> None:
    assert _outcome(split_items, raw, quotes=quotes) == _outcome(_ref_split_items, raw, quotes=quotes)
