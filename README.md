# pefftacular

[![PyPI](https://img.shields.io/pypi/v/pefftacular)](https://pypi.org/project/pefftacular/)
[![Python Package](https://github.com/tacular-omics/pefftacular/actions/workflows/python-package.yml/badge.svg)](https://github.com/tacular-omics/pefftacular/actions/workflows/python-package.yml)
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

**read_peff** — load everything into memory at once:

```python
from pefftacular import read_peff

header, entries = read_peff("proteins.peff")

for entry in entries:
    print(entry.db_unique_id, entry.pname, len(entry.sequence))
```

**PeffReader** — iterate lazily without loading the full file:

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
| `variant_simple` | `tuple[VariantSimple, ...]` | Simple sequence variants |
| `variant_complex` | `tuple[VariantComplex, ...]` | Multi-residue variants (start, end, new sequence, optional tag) |
| `mod_res_unimod` | `tuple[ModResUnimod, ...]` | UniMod modification sites |
| `mod_res_psi` | `tuple[ModResPsi, ...]` | PSI-MOD modification sites |
| `mod_res` | `tuple[ModRes, ...]` | Other named modification sites |
| `processed` | `tuple[Processed, ...]` | Processed sequence forms |
| `custom_values` | `dict[str, tuple[CustomKeyValue, ...]]` | Header-declared custom keys, parsed by their `CustomKeyDef` |
| `extra` | `dict[str, str]` | Non-standard keys with no `CustomKeyDef` |

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

**Custom keys (declared via `# CustomKeyDef=` in the header):**

When the database header declares a custom key, entry values for that key are
parsed using its `RegExp` / `FieldNames` / `FieldTypes` and exposed as typed
fields on `entry.custom_values`. The original item text is preserved in `raw`
for lossless round-trips.

Header excerpt:

```
# CustomKeyDef=(KeyName=SecondaryStructure|Description="..."|ConceptCURIE=BAO:0000014|RegExp="([0-9]+)\|([0-9]+)\|([A-Za-z]+:[0-9]+)?\|(.+)"|FieldNames=StartPosition,EndPosition,CURIE,Description|FieldTypes=integer,integer,string,string)
```

Entry usage:

```
>cu:P00001 \SecondaryStructure=(10|20|ncithesaurus:C47937|Helix)
```

Access:

```python
ss = entry.custom_values["SecondaryStructure"]
ss[0].fields["StartPosition"]    # 10 (int)
ss[0].fields["Description"]      # "Helix"
```

Supported `FieldTypes` are XSD basic types (`string`, `integer`, `decimal`,
`boolean`, `date`, `time`) plus `enumeration(a|b|c)`. Coercion failures and
enumeration mismatches emit `UserWarning` and fall back to the raw string.
If no `RegExp` is declared, the value is split on `|` and zipped with
`FieldNames`.

**Other non-standard keys** (no `CustomKeyDef` registered) still land in
`entry.extra` as raw strings:

```python
value = entry.extra.get("MyCustomKey")
```

## Writing

Build a header and entries, then write:

```python
from pefftacular import DatabaseHeader, FileHeader, SequenceEntry, write_peff

db_header = DatabaseHeader(
    prefix="sp",
    db_name="SwissProt",
    db_version="2024_01",
    number_of_entries=1,
)

file_header = FileHeader(
    peff_version="1.0",
    databases=(db_header,),
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

`dest` can be a file path string, a `pathlib.Path`, or a text-mode file object.

## Error handling

Parse errors raise `PeffParseError`:

```python
from pefftacular import PeffParseError, read_peff

try:
    header, entries = read_peff("malformed.peff")
except PeffParseError as e:
    print(e.line)     # the offending line number
    print(e.context)  # surrounding context string
```

Write errors raise `PeffWriteError`:

```python
from pefftacular import PeffWriteError

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
