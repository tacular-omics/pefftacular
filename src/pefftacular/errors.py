"""PEFF-specific error types."""


class PeffParseError(ValueError):
    """Raised when PEFF input cannot be parsed."""

    def __init__(self, message: str, *, line: int | None = None, context: str | None = None) -> None:
        self.line = line
        self.context = context
        super().__init__(message if line is None else f"Line {line}: {message}")


class PeffWriteError(ValueError):
    """Raised when a model object cannot be serialized to PEFF."""
