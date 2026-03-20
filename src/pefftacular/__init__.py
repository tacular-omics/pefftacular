"""pefftacular — A pure-Python PEFF parsing and writing library."""

from pefftacular._models import (
    CustomKeyDef,
    DatabaseHeader,
    DisulfideBond,
    FileHeader,
    ModRes,
    ModResPsi,
    ModResUnimod,
    OptionalTagDef,
    Processed,
    Proteoform,
    SequenceEntry,
    SequenceRange,
    VariantComplex,
    VariantSimple,
)
from pefftacular._parser import PeffReader, read_peff
from pefftacular._writer import write_peff
from pefftacular.errors import PeffParseError, PeffWriteError

__all__ = [
    "CustomKeyDef",
    "DatabaseHeader",
    "DisulfideBond",
    "FileHeader",
    "ModRes",
    "ModResPsi",
    "ModResUnimod",
    "OptionalTagDef",
    "PeffParseError",
    "PeffReader",
    "PeffWriteError",
    "Processed",
    "Proteoform",
    "SequenceEntry",
    "SequenceRange",
    "VariantComplex",
    "VariantSimple",
    "read_peff",
    "write_peff",
]
