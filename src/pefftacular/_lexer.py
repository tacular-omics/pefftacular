"""Depth-tracking tokenizers for PEFF description lines."""

from pefftacular.errors import PeffParseError


def split_items(raw: str) -> list[str]:
    """Split parenthesized multi-item values into individual items.

    ``(A|B)(C|D)`` -> ``["A|B", "C|D"]``
    Single unparenthesized values are returned as-is: ``"110"`` -> ``["110"]``

    Raises:
        PeffParseError: On mismatched parentheses.
    """
    if not raw:
        return [raw]
    if not raw.startswith("("):
        if ")" in raw:
            raise PeffParseError("Unexpected ')' in value", context=raw)
        return [raw]

    items: list[str] = []
    depth = 0
    item_start = -1

    for i, ch in enumerate(raw):
        match ch:
            case "(":
                if depth == 0:
                    item_start = i + 1
                depth += 1
            case ")":
                depth -= 1
                if depth == 0:
                    items.append(raw[item_start:i])
                elif depth < 0:
                    raise PeffParseError(f"Unexpected ')' at position {i}", context=raw)

    if depth != 0:
        raise PeffParseError("Unclosed '(' in value", context=raw)

    return items


def split_fields(item: str) -> list[str]:
    """Split on ``|`` at paren depth 0 only, preserving empty strings."""
    fields: list[str] = []
    depth = 0
    start = 0

    for i, ch in enumerate(item):
        match ch:
            case "(":
                depth += 1
            case ")":
                depth -= 1
            case "|" if depth == 0:
                fields.append(item[start:i])
                start = i + 1

    fields.append(item[start:])
    return fields


def split_description_keys(rest: str) -> dict[str, str]:
    r"""Parse the annotation portion of a description line into ``{key: raw_value}``.

    Scans for ``\Key=value`` boundaries at paren depth 0.
    A new key begins when ``\`` is encountered at depth 0 and is either at
    the start of the string or preceded by a space.
    """
    if not rest:
        return {}

    keys: dict[str, str] = {}
    depth = 0
    current_start: int | None = None  # index of the '\' that starts the current token

    for i, ch in enumerate(rest):
        match ch:
            case "(":
                depth += 1
            case ")":
                depth -= 1
            case "\\" if depth == 0 and (i == 0 or rest[i - 1] == " "):
                # Close previous key-value if one is open
                if current_start is not None:
                    _store_key_value(keys, rest[current_start:i].rstrip())
                current_start = i

    # Store last token
    if current_start is not None:
        _store_key_value(keys, rest[current_start:])

    return keys


def _store_key_value(keys: dict[str, str], token: str) -> None:
    r"""Parse a ``\Key=value`` token and insert into *keys*."""
    # Strip leading backslash
    token = token.lstrip("\\")
    eq_idx = token.find("=")
    if eq_idx == -1:
        # Key with no value (e.g. ``\Decoy``)
        keys[token] = ""
    else:
        keys[token[:eq_idx]] = token[eq_idx + 1 :]
