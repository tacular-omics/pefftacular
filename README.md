# pefftacular

[![PyPI](https://img.shields.io/pypi/v/pefftacular)](https://pypi.org/project/pefftacular/)
[![CI](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/tacular-omics/pefftacular)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/pefftacular)](https://pypi.org/project/pefftacular/)

Python library for reading and writing [PEFF](http://www.psidev.info/peff) (PSI Extended FASTA Format) files. PEFF is a superset of FASTA used in proteomics that carries rich per-entry annotations — PTMs, variants, processed forms, and more — encoded directly in the sequence header.

## Install

```bash
pip install pefftacular
```

Dev install:

```bash
just install
```

## Quick start

**Small files** — load everything into memory at once:

```python
from pefftacular import read_peff

header, entries = read_peff("proteins.peff")

for entry in entries:
    print(entry.db_unique_id, entry.pname, len(entry.sequence))
```

**Large files** — iterate lazily without loading the full file:

```python
from pefftacular import PeffReader

with PeffReader("proteins.peff") as reader:
    file_header = reader.header
    for entry in reader:
        process(entry)
```

## Data model

`read_peff` and `PeffReader` yield `SequenceEntry` objects with these fields:

| Field | Type | Description |
|---|---|---|
| `prefix` | `str` | Database prefix (e.g. `sp`, `tr`) |
| `db_unique_id` | `str` | Accession (e.g. `P12345`) |
| `sequence` | `str` | Amino acid sequence |
| `pname` | `str \| None` | Protein name (`\\PName=`) |
| `gname` | `str \| None` | Gene name (`\\GName=`) |
| `ncbi_tax_id` | `int \| None` | NCBI taxonomy ID (`\\NcbiTaxId=`) |
| `length` | `int \| None` | Sequence length (`\\Length=`) |
| `sv` | `int \| None` | Sequence version (`\\SV=`) |
| `ev` | `int \| None` | Entry version (`\\EV=`) |
| `pe` | `int \| None` | Protein existence level (`\\PE=`) |
| `variant_simple` | `list[VariantSimple]` | Simple sequence variants |
| `variant_complex` | `list[str]` | Complex variant strings (unparsed) |
| `mod_res_unimod` | `list[ModResUnimod]` | UniMod modification sites |
| `mod_res_psi` | `list[ModResPsi]` | PSI-MOD modification sites |
| `mod_res` | `list[ModRes]` | Other named modification sites |
| `processed` | `list[Processed]` | Processed sequence forms |
| `extra` | `dict[str, str]` | Non-standard key/value pairs |

## Annotations

**Variants:**

```python
from pefftacular import read_peff

_, entries = read_peff("proteins.peff")
entry = entries[0]

for v in entry.variant_simple:
    print(v.position, v.new_amino_acid, v.tag)
    # e.g. 42, "K", "rs12345"
```

**Modifications (UniMod):**

```python
for mod in entry.mod_res_unimod:
    print(mod.position, mod.accession, mod.name)
    # e.g. 17, "21", "Phospho"
```

**Modifications (PSI-MOD):**

```python
for mod in entry.mod_res_psi:
    print(mod.position, mod.accession, mod.name)
    # e.g. 17, "MOD:00696", "phosphorylated residue"
```

**Processed forms:**

```python
for proc in entry.processed:
    print(proc.start_pos, proc.end_pos, proc.accession, proc.name)
    # e.g. 1, 24, "PRO_0000012345", "Signal peptide"
```

**Non-standard keys:**

```python
value = entry.extra.get("MyCustomKey")
```

## Writing

Build a header and entries, then write:

```python
from pefftacular import write_peff
from pefftacular.models import FileHeader, DatabaseHeader, SequenceEntry

db_header = DatabaseHeader(
    prefix="sp",
    db_name="SwissProt",
    db_version="2024_01",
    number_of_entries=1,
)

file_header = FileHeader(
    version="1.0",
    databases=[db_header],
)

entry = SequenceEntry(
    prefix="sp",
    db_unique_id="P12345",
    sequence="MKTIIALSYIFCLVFA",
    pname="Example protein",
    gname="EXMP",
)

write_peff(file_header, [entry], "output.peff")
```

`dest` can be a file path string or a `pathlib.Path`. Pass an open binary file object to write to an existing stream.

## Error handling

Parse errors raise `PeffParseError`:

```python
from pefftacular import read_peff
from pefftacular.exceptions import PeffParseError

try:
    header, entries = read_peff("malformed.peff")
except PeffParseError as e:
    print(e.line)     # the offending line number
    print(e.context)  # surrounding context string
```

Write errors raise `PeffWriteError`:

```python
from pefftacular.exceptions import PeffWriteError

try:
    write_peff(file_header, entries, "/read-only/output.peff")
except PeffWriteError as e:
    print(e)
```

## Development

```bash
just install      # install dependencies
just test         # run tests
just test-v       # run tests (verbose)
just test-file tests/test_reader.py   # run a single test file
just cov          # run tests with coverage
just lint         # ruff lint
just format       # ruff format
just check        # lint + type check + test
just build        # build the package
just clean        # remove cache files
just docs         # serve docs locally
just docs-deploy  # deploy docs to GitHub Pages
```

## License

[MIT](LICENSE)
