"""Depth-tracking tokenizers for PEFF description lines."""

from pefftacular.errors import PeffParseError


def _unescape_component(s: str) -> str:
    r"""Reverse PEFF backslash-escaping in a single component.

    Per the spec (section 3.3.3), a ``\`` escapes a following ``|``, ``(``,
    ``)`` or ``\`` so the character is taken literally rather than as a
    separator / paren. Backslashes before any other character are preserved.
    """
    if "\\" not in s:
        return s
    out: list[str] = []
    i = 0
    length = len(s)
    while i < length:
        ch = s[i]
        if ch == "\\" and i + 1 < length and s[i + 1] in ("|", "(", ")", "\\"):
            out.append(s[i + 1])
            i += 2
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def _has_unescaped(s: str, target: str) -> bool:
    """Whether *target* appears in *s* not preceded by an escaping backslash."""
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


def split_items(raw: str) -> list[str]:
    r"""Split parenthesized multi-item values into individual items.

    ``(A|B)(C|D)`` -> ``["A|B", "C|D"]``
    Single unparenthesized values are returned as-is: ``"110"`` -> ``["110"]``

    Backslash-escaped parens (``\(`` / ``\)``) do not affect nesting depth, so
    an item may contain an escaped unpaired paren. The returned substrings keep
    their escapes; unescaping happens per-component in :func:`split_fields`.

    Raises:
        PeffParseError: On mismatched parentheses.
    """
    if not raw:
        return [raw]
    if not raw.startswith("("):
        if _has_unescaped(raw, ")"):
            raise PeffParseError(
                "Unexpected ')' in value",
                context=raw,
                hint=r"Escape a literal ')' as '\)', or wrap multi-item values in matching parentheses",
            )
        return [raw]

    items: list[str] = []
    depth = 0
    item_start = -1
    i = 0
    length = len(raw)

    while i < length:
        ch = raw[i]
        if ch == "\\" and i + 1 < length:
            i += 2
            continue
        if ch == "(":
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


def split_fields(item: str, *, unescape: bool = False) -> list[str]:
    r"""Split an item on ``|`` at paren depth 0 into its components.

    Two escaping regimes are supported:

    * Default (``unescape=False``) — used for ``CustomKeyDef`` header values,
      where ``|`` inside ``"..."`` quoted spans is ignored and backslashes are
      preserved verbatim (so embedded regexes keep their meaning).
    * ``unescape=True`` — used for entry description items, where the spec's
      backslash escaping applies: ``\|`` is a literal pipe (not a separator),
      ``\(`` / ``\)`` are literal parens (no depth change), and each returned
      component is unescaped via :func:`_unescape_component`. Quotes carry no
      special meaning in this mode.

    Empty strings between separators are preserved in both modes.
    """
    if unescape:
        return _split_fields_escaped(item)

    fields: list[str] = []
    depth = 0
    start = 0
    in_quote = False
    i = 0
    length = len(item)

    while i < length:
        ch = item[i]
        if in_quote:
            if ch == "\\" and i + 1 < length:
                i += 2
                continue
            if ch == '"':
                in_quote = False
        else:
            if ch == '"':
                in_quote = True
            elif ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "|" and depth == 0:
                fields.append(item[start:i])
                start = i + 1
        i += 1

    fields.append(item[start:])
    return fields


def _split_fields_escaped(item: str) -> list[str]:
    r"""Split entry-item components on unescaped ``|``, then unescape each.

    Backslash escapes suppress separator/paren meaning: ``\|`` stays within a
    component and ``\(`` / ``\)`` do not change nesting depth.
    """
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
            fields.append(_unescape_component(item[start:i]))
            start = i + 1
        i += 1

    fields.append(_unescape_component(item[start:]))
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
