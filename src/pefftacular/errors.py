"""PEFF-specific error and warning types.

Every exception raised by pefftacular derives from :class:`PeffError`, and every
non-fatal spec violation is emitted through :class:`PeffWarning`. This lets
callers — and automated coding agents — catch or filter *all* library-specific
signals with a single clause::

    try:
        header, entries = read_peff(path)
    except PeffError as err:
        ...  # any pefftacular parse/write failure

    warnings.simplefilter("error", PeffWarning)  # turn spec violations into errors

Parse failures carry machine-readable attributes (``line``, ``context``,
``hint``) in addition to a human-readable message, and attach the offending text
and a remediation ``hint`` as exception *notes* so they show up in tracebacks.
"""

from __future__ import annotations


class PeffError(ValueError):
    """Base class for all pefftacular errors.

    Subclasses :class:`ValueError` for backwards compatibility with callers that
    catch ``ValueError``.
    """


class PeffParseError(PeffError):
    """Raised when PEFF input cannot be parsed.

    Attributes:
        line: 1-based absolute line number the failure was detected on, if known.
        context: The offending text (e.g. the raw item that failed to parse).
        hint: A short, actionable suggestion for how to fix the input.
    """

    def __init__(
        self,
        message: str,
        *,
        line: int | None = None,
        context: str | None = None,
        hint: str | None = None,
    ) -> None:
        self.line = line
        self.context = context
        self.hint = hint
        super().__init__(message if line is None else f"Line {line}: {message}")
        # Notes surface in tracebacks (PEP 678) without changing ``str(err)``,
        # so agents reading a stack trace see the offending text and the fix.
        if context is not None:
            self.add_note(f"offending text: {context!r}")
        if hint is not None:
            self.add_note(f"hint: {hint}")


class PeffWriteError(PeffError):
    """Raised when a model object cannot be serialized to PEFF.

    Attributes:
        hint: A short, actionable suggestion for how to fix the model.
    """

    def __init__(self, message: str, *, hint: str | None = None) -> None:
        self.hint = hint
        super().__init__(message)
        if hint is not None:
            self.add_note(f"hint: {hint}")


class PeffWarning(UserWarning):
    """Category for non-fatal PEFF spec violations detected while reading.

    Parsing stays permissive: the data is always returned, but anything that
    violates a spec ``MUST`` rule (out-of-range positions, missing required
    fields, header/entry-count mismatches, un-coercible custom values, …) is
    reported through this category. Because it subclasses :class:`UserWarning`,
    existing ``UserWarning`` filters keep working; agents that want to treat PEFF
    problems as hard errors can ``warnings.simplefilter("error", PeffWarning)``.
    """
