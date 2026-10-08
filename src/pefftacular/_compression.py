"""The ``compression`` argument shared by the readers and the writer."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Literal, cast, get_args

from pefftacular.errors import PeffError

Compression = Literal["infer", "gzip", "bz2", "xz"] | None
"""How a PEFF file is (de)compressed, like the ``compression`` argument of pandas.

- ``"infer"`` (the default): readers detect gzip, bzip2 or xz from the magic bytes;
  writers compress by the path suffix (``.gz``, ``.bz2``, ``.xz``, any case). An open
  handle is read and written as plain text.
- ``"gzip"``, ``"bz2"``, ``"xz"``: always that format, whatever the suffix or content.
- ``None``: never compressed.
"""

_VALID: tuple[str | None, ...] = get_args(Literal["infer", "gzip", "bz2", "xz"]) + (None,)
_SUFFIXES = {".gz": "gzip", ".bz2": "bz2", ".xz": "xz"}
# Leading bytes of each format.
MAGIC = {"gzip": b"\x1f\x8b", "bz2": b"BZh", "xz": b"\xfd7zXZ\x00"}
MODULES = {"gzip": "gzip", "bz2": "bz2", "xz": "lzma"}


def check_compression(compression: object) -> Compression:
    """Return ``compression`` if it is a valid value; raise ``PeffError`` otherwise."""
    if compression is None or (isinstance(compression, str) and compression in _VALID):
        return cast("Compression", compression)
    err = PeffError(f"Unknown compression {compression!r}; valid values are {', '.join(map(repr, _VALID))}")
    err.add_note('hint: use "infer" (the default), "gzip", "bz2", "xz" or None')
    raise err


def from_suffix(path: Path) -> str | None:
    """The compression a path's suffix names (``.gz``, ``.bz2``, ``.xz``, any case), else ``None``."""
    return _SUFFIXES.get(path.suffix.lower())


def sniff(head: bytes) -> str | None:
    """The compression whose magic bytes ``head`` starts with, else ``None``."""
    return next((kind for kind, magic in MAGIC.items() if head.startswith(magic)), None)


def is_binary(handle: object) -> bool:
    """Whether an open handle reads/writes ``bytes`` rather than ``str``."""
    if isinstance(handle, io.TextIOBase):
        return False
    if isinstance(handle, (io.RawIOBase, io.BufferedIOBase)):
        return True
    mode = getattr(handle, "mode", None)
    return isinstance(mode, str) and "b" in mode
