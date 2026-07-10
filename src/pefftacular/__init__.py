"""pefftacular — A pure-Python PEFF parsing and writing library."""

import logging

from pefftacular._models import (
    CustomKeyDef,
    CustomKeyValue,
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
from pefftacular.errors import (
    PeffError,
    PeffParseError,
    PeffWarning,
    PeffWriteError,
)

# Library logging follows the stdlib convention: attach a NullHandler so importing
# pefftacular never emits output on its own. Applications (and coding agents that
# want a behavioural trace) opt in with, e.g.::
#
#     import logging
#     logging.basicConfig(level=logging.DEBUG)
#     logging.getLogger("pefftacular").setLevel(logging.DEBUG)
logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "CustomKeyDef",
    "CustomKeyValue",
    "DatabaseHeader",
    "DisulfideBond",
    "FileHeader",
    "ModRes",
    "ModResPsi",
    "ModResUnimod",
    "OptionalTagDef",
    "PeffError",
    "PeffParseError",
    "PeffReader",
    "PeffWarning",
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
