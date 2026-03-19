"""pefftacular — A pure-Python PEFF parsing and writing library."""

from pefftacular._models import (
    CustomKeyDef,
    DatabaseHeader,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    Processed,
    SequenceEntry,
    VariantComplex,
    VariantSimple,
)
from pefftacular._parser import PeffReader, read_peff
from pefftacular._writer import write_peff
from pefftacular.errors import PeffParseError, PeffWriteError

__all__ = [
    "CustomKeyDef",
    "DatabaseHeader",
    "FileHeader",
    "ModRes",
    "ModResPsi",
    "ModResUnimod",
    "PeffParseError",
    "PeffReader",
    "PeffWriteError",
    "Processed",
    "SequenceEntry",
    "VariantComplex",
    "VariantSimple",
    "read_peff",
    "write_peff",
]
