# Getting started

## Installation

pefftacular needs Python 3.12 or later and has no runtime dependencies.

```bash
pip install pefftacular
# or, in a uv project
uv add pefftacular
```

## A file to work with

The examples on this page use this small PEFF file. Save it as `proteins.peff`:

```text
# PEFF 1.0
# GeneralComment=Two-entry example for the pefftacular docs
# //
# DbName=SwissProt
# Prefix=sp
# DbVersion=2026_03
# DbSource=https://www.uniprot.org
# NumberOfEntries=2
# SequenceType=AA
# //
>sp:P69905 \PName=Hemoglobin subunit alpha \GName=HBA1 \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=30 \SV=2 \PE=1 \ModResPsi=(4|MOD:00046|O-phospho-L-serine)
MVLSPADKTNVKAAWGKVGAHAGEYGAEAL
>sp:P68871 \PName=Hemoglobin subunit beta \GName=HBB \NcbiTaxId=9606 \TaxName=Homo sapiens \Length=30 \SV=2 \PE=1 \VariantSimple=(7|V|sickle cell)
MVHLTPEEKSAVTALWGKVNVDEVGGEALG
```

## Reading

### Everything at once: `read_peff`

`read_peff` returns the file header and a list of entries:

```python
from pefftacular import read_peff

header, entries = read_peff("proteins.peff")

print(len(entries))
# 2
for entry in entries:
    print(entry.prefix, entry.db_unique_id, entry.gname, entry.pname)
# sp P69905 HBA1 Hemoglobin subunit alpha
# sp P68871 HBB Hemoglobin subunit beta
```

### One entry at a time: `PeffReader`

`PeffReader` parses the header up front and then yields entries lazily, so memory use stays flat
on large databases. It must be used as a context manager: the file is opened on entering the
`with` block and closed on leaving it, and using the reader outside one raises `RuntimeError`:

```python
from pefftacular import PeffReader

with PeffReader("proteins.peff") as reader:
    print(reader.header.databases[0].db_name)
    for entry in reader:
        print(entry.db_unique_id, len(entry.sequence))
# SwissProt
# P69905 30
# P68871 30
```

### Paths, file objects and strings

Both readers accept a path (`str` or `pathlib.Path`) or any text-mode file object. A `str` is
always treated as a **path**. To parse PEFF text you already hold in memory, wrap it in
`io.StringIO`:

```python
import io

from pefftacular import read_peff

text = """# PEFF 1.0
# //
# DbName=Demo
# Prefix=db
# DbVersion=1
# DbSource=local
# NumberOfEntries=1
# SequenceType=AA
# //
>db:X1 \\PName=Demo protein
PEPTIDE
"""

header, entries = read_peff(io.StringIO(text))
print(entries[0].pname, entries[0].sequence)
# Demo protein PEPTIDE
```

## The data model

Everything pefftacular returns is a frozen dataclass: immutable and comparable with `==`.
`SequenceEntry` is **not hashable**, because `custom_values` and `extra` are dicts; use
`(entry.prefix, entry.db_unique_id)` as a key instead.

```text
FileHeader
├── peff_version          "1.0"
├── general_comments      file-level "# GeneralComment=" lines
└── databases             one DatabaseHeader per "# //" block
    └── DatabaseHeader
        ├── prefix, db_name, db_version, db_sources, number_of_entries, ...
        ├── custom_key_defs   CustomKeyDef, one per "# CustomKeyDef=" line
        └── optional_tag_defs OptionalTagDef, one per "# OptionalTagDef=" line

SequenceEntry             one per ">" description line
├── prefix, db_unique_id, sequence
├── pname, gname, ncbi_tax_id, tax_name, length, sv, ev, pe, ...
├── variant_simple, variant_complex, mod_res_unimod, mod_res_psi, mod_res,
│   processed, disulfide_bond, proteoform      (tuples of annotation objects)
├── custom_values         keys declared by a CustomKeyDef
└── extra                 any other \Key=value pair, as raw strings
```

### The file header

```python
from pefftacular import read_peff

header, entries = read_peff("proteins.peff")

print(header.peff_version)
# 1.0
print(header.general_comments)
# ('Two-entry example for the pefftacular docs',)

db = header.databases[0]
print(db.prefix, db.db_name, db.db_version, db.number_of_entries)
# sp SwissProt 2026_03 2
print(db.db_sources, db.sequence_type)
# ('https://www.uniprot.org',) AA
```

A PEFF file can hold several databases, each with its own prefix. Every entry's `prefix` names
the database it belongs to, so you can match them up:

```python
dbs = {db.prefix: db for db in header.databases}
for entry in entries:
    print(entry.db_unique_id, "from", dbs[entry.prefix].db_name)
# P69905 from SwissProt
# P68871 from SwissProt
```

Header keys pefftacular does not model (for example `SpecificKey` blocks) are kept as raw
strings in `DatabaseHeader.extra`.

### Sequence entries

The description line `>sp:P69905 \PName=... \GName=HBA1 ...` becomes a `SequenceEntry`.
The part before the colon is `prefix`, the part after is `db_unique_id`, and each standard key
maps to a typed field:

| PEFF key | Field | Type |
|---|---|---|
| (prefix) | `prefix` | `str` |
| (unique id) | `db_unique_id` | `str` |
| (sequence lines) | `sequence` | `str` |
| `\ID` | `id` | `str \| None` |
| `\DbUniqueId` | `db_unique_id_key` | `str \| None` |
| `\PName` | `pname` | `str \| None` |
| `\GName` | `gname` | `str \| None` |
| `\NcbiTaxId` (or `\OX`) | `ncbi_tax_id` | `int \| None` |
| `\TaxName` | `tax_name` | `str \| None` |
| `\Length` | `length` | `int \| None` |
| `\SV` / `\EV` / `\PE` | `sv` / `ev` / `pe` | `int \| None` |
| `\Decoy` | `decoy` | `bool \| None` |
| `\Comment` | `comment` | `str \| None` |
| `\VariantSimple` | `variant_simple` | `tuple[VariantSimple, ...]` |
| `\VariantComplex` | `variant_complex` | `tuple[VariantComplex, ...]` |
| `\ModResUnimod` | `mod_res_unimod` | `tuple[ModResUnimod, ...]` |
| `\ModResPsi` | `mod_res_psi` | `tuple[ModResPsi, ...]` |
| `\ModRes` | `mod_res` | `tuple[ModRes, ...]` |
| `\Processed` | `processed` | `tuple[Processed, ...]` |
| `\DisulfideBond` | `disulfide_bond` | `tuple[DisulfideBond, ...]` |
| `\Proteoform` | `proteoform` | `tuple[Proteoform, ...]` |
| custom keys | `custom_values` | `dict[str, tuple[CustomKeyValue, ...]]` |
| anything else | `extra` | `dict[str, str]` |

```python
alpha = entries[0]
print(alpha.ncbi_tax_id, alpha.tax_name, alpha.length, alpha.sv, alpha.pe)
# 9606 Homo sapiens 30 2 1
print(alpha.mod_res_psi[0])
# ModResPsi(positions=(4,), accession='MOD:00046', name='O-phospho-L-serine', tag=None, annot_id=None)
```

Annotation types are covered in detail in the [annotations guide](annotations.md).

### Changing an entry

Models are frozen, so "editing" means making a copy with `dataclasses.replace`:

```python
from dataclasses import replace

renamed = replace(alpha, gname="HBA2")
print(alpha.gname, renamed.gname)
# HBA1 HBA2
```

## Writing

`write_peff(header, entries, dest)` writes a complete file. `dest` is a path or a text-mode file
object; `entries` can be any iterable, including a generator.

```python
from pefftacular import DatabaseHeader, FileHeader, SequenceEntry, write_peff

header = FileHeader(
    peff_version="1.0",
    databases=(
        DatabaseHeader(
            prefix="my",
            db_name="MyProteins",
            db_version="1.0",
            db_sources=("in-house",),
            number_of_entries=1,
            sequence_type="AA",
        ),
    ),
)

entry = SequenceEntry(
    prefix="my",
    db_unique_id="PROT001",
    sequence="MKTIIALSYIFCLVFA",
    pname="Example protein",
    gname="EXMP",
    length=16,
)

write_peff(header, [entry], "output.peff")
print(open("output.peff").read())
```

```text
# PEFF 1.0
# //
# DbName=MyProteins
# Prefix=my
# DbVersion=1.0
# DbSource=in-house
# NumberOfEntries=1
# SequenceType=AA
# //
>my:PROT001 \Length=16 \PName=Example protein \GName=EXMP
MKTIIALSYIFCLVFA
```

The writer:

- emits keys in a fixed canonical order, so output is stable across runs,
- wraps sequences at 60 residues per line,
- backslash-escapes `\`, `|` and unbalanced parentheses inside annotation fields.

It does **not** fill in `Length` or `NumberOfEntries` for you. Set them yourself if you want
them in the file; the reader warns if they disagree with the data.

### Filtering a file

Read, transform, write. `PeffReader` streams the input; `write_peff` collects the entries it is
given into a list before writing, so the entries you keep must fit in memory:

```python
from pefftacular import PeffReader, write_peff

with PeffReader("proteins.peff") as reader:
    header = reader.header
    kept = [e for e in reader if e.variant_simple]

write_peff(header, kept, "with_variants.peff")
print([e.db_unique_id for e in kept])
# ['P68871']
```

If you filter entries out, the header's `NumberOfEntries` still says the old count and the
reader will warn when it reads the new file. Update it with `dataclasses.replace`:

```python
from dataclasses import replace

db = replace(header.databases[0], number_of_entries=len(kept))
write_peff(replace(header, databases=(db,)), kept, "with_variants.peff")
```

## Errors and warnings

### Parse errors

Input that cannot be parsed at all raises `PeffParseError`. It carries `line`, `context` (the
offending text) and `hint`. Every pefftacular exception derives from `PeffError`, which is a
`ValueError`.

```python
import io

from pefftacular import PeffError, PeffParseError, read_peff

try:
    read_peff(io.StringIO("#PEFF 1.0\n>sp:P1\nMK\n"))
except PeffParseError as err:
    print(err.line)
    print(err.context)
    print(err.hint)
# 1
# #PEFF 1.0
# The first non-blank line must be exactly '# PEFF 1.0' (note the space after '#')
```

The context and hint are also attached as exception notes, so they appear in a traceback.
`write_peff` raises `PeffWriteError` for entries it cannot serialize, such as an empty sequence.

### Spec-violation warnings

Reading is permissive. If a file breaks a PEFF `MUST` rule but can still be parsed (a missing
mandatory header key, a `\Length` that does not match the sequence, an annotation position
outside the sequence, a wrong `NumberOfEntries`), you get the data plus a `PeffWarning`:

```python
import io
import warnings

from pefftacular import PeffWarning, read_peff

bad = """# PEFF 1.0
# //
# DbName=Demo
# Prefix=db
# DbVersion=1
# DbSource=local
# NumberOfEntries=1
# SequenceType=AA
# //
>db:X1 \\Length=99
PEPTIDE
"""

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    header, entries = read_peff(io.StringIO(bad))

print(entries[0].sequence)
# PEPTIDE
for w in caught:
    print(w.category.__name__, w.message)
# PeffWarning Entry 'X1': Length=99 but sequence has 7 residues
```

For strict parsing, turn the warnings into errors:

```python
warnings.simplefilter("error", PeffWarning)
try:
    read_peff(io.StringIO(bad))
except PeffWarning as w:
    print("rejected:", w)
# rejected: Entry 'X1': Length=99 but sequence has 7 residues
```

`PeffWarning` subclasses `UserWarning`, so existing `UserWarning` filters still apply. The
deprecated `\Variant=` key also warns with `PeffWarning`, and its value lands in `entry.extra`.

## Logging

pefftacular never configures logging itself. To see what it is doing, enable the
`pefftacular` logger:

```python
import logging

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("pefftacular").setLevel(logging.DEBUG)
```

File opens, header parsing and entry counts log at `DEBUG`; totals read and written log at
`INFO`, under the `pefftacular.parser` and `pefftacular.writer` loggers.
