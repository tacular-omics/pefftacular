# pefftacular

[![PyPI](https://img.shields.io/pypi/v/pefftacular)](https://pypi.org/project/pefftacular/)
[![Python Package](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/tacular-omics/pefftacular)](https://github.com/tacular-omics/pefftacular/blob/main/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/pefftacular)](https://pypi.org/project/pefftacular/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22925639.svg)](https://doi.org/10.5281/zenodo.22925639)
[![Docs](https://img.shields.io/badge/docs-tacular--omics.github.io-blue)](https://tacular-omics.github.io/pefftacular/)

A pure-Python library for reading and writing [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format) files — the proteomics community's FASTA extension for carrying rich per-entry annotations (PTMs, sequence variants, processed forms, and more) directly in the sequence header. It's for anyone building or consuming protein sequence databases that need more than a bare FASTA header can hold.

pefftacular parses PEFF into typed, structured objects instead of leaving you to regex the header yourself, and it's permissive by default: malformed-but-recoverable files still parse, with spec violations reported as warnings rather than hard failures.

**Documentation:** <https://tacular-omics.github.io/pefftacular/>

## Highlights

- **Zero dependencies** — pure Python, nothing else to install.
- **Two ways to read** — `read_peff` for the whole file at once, `PeffReader` to stream entries lazily without loading everything into memory.
- **Rich annotations as typed data** — variants, UniMod/PSI-MOD modification sites, processed forms, and header-declared custom keys all come back as structured fields, not strings you parse yourself.
- **Permissive reading, strict when you want it** — spec violations are reported through `PeffWarning` (opt in to `warnings.simplefilter("error", PeffWarning)` for strict parsing) while `PeffParseError` still carries a line, context, and a repair hint for structural failures.
- **Round-trip safe writing**, including the escaping rules PEFF requires in description lines.
- **Shares its API shape with [fastatacular](https://github.com/tacular-omics/fastatacular)**, the plain-FASTA sibling library, so switching formats doesn't mean relearning the interface.

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

Both accept gzip, bzip2 and xz compressed files (`proteins.peff.gz`) directly.

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
    print(mod.positions, mod.accession, mod.name)
    # e.g. (17,), "21", "Phospho"
```

**Modifications (PSI-MOD):**

```python
for mod in entry.mod_res_psi:
    print(mod.positions, mod.accession, mod.name)
    # e.g. (17,), "MOD:00696", "phosphorylated residue"
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
enumeration mismatches emit `PeffWarning` (a `UserWarning` subclass) and fall back to the raw string.
If no `RegExp` is declared, the value is split on `|` and zipped with
`FieldNames`.

**Other non-standard keys** (no `CustomKeyDef` registered) still land in
`entry.extra` as raw strings:

```python
value = entry.extra.get("MyCustomKey")
```

## FASTA and ProForma

```python
from pefftacular import SequenceEntry

entry = SequenceEntry.from_fasta(
    "sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens OX=9606 GN=YWHAB PE=1 SV=3",
    "MTMDKSELVQKAKLAEQAERYDDMAAAMK",
)
entry.prefix, entry.db_unique_id, entry.id   # ("sp", "P31946", "1433B_HUMAN")
header, sequence = entry.to_fasta()          # back to a plain FASTA record

entry.to_proforma()                  # \ModResPsi sites, e.g. "MS[MOD:00046]TK..."
entry.to_proforma(mods="unimod")     # from \ModResUnimod sites
entry.to_proforma(variants=entry.variant_simple[:1])   # with a substitution applied
```

These use plain strings (pefftacular has no dependencies). With fastatacular, pass
`(record.raw_header, record.sequence)`. See [From FASTA](docs/fasta.md) and
[ProForma strings](docs/annotations.md#proforma-strings).

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
`entries` can be any iterable, including a generator: it is consumed once, as a stream.

Every written line is parsed back and compared with the entry, so a value that would
read back differently raises `PeffWriteError` instead of corrupting the file. For
entries that came unchanged from `read_peff`, pass `verify=False` to skip that check
(about 3x faster writes):

```python
header, entries = read_peff("proteins.peff.gz")
write_peff(header, (e for e in entries if e.variant_simple), "variants.peff", verify=False)
```

## Error handling

Every exception derives from `PeffError` (a `ValueError` subclass), so you can
catch any failure with one clause. Parse errors carry structured, actionable
detail — `.line`, `.context`, and a `.hint` — and attach the offending text and
the hint as exception *notes*, so they also show up in tracebacks:

```python
from pefftacular import PeffError, PeffParseError, read_peff

try:
    header, entries = read_peff("malformed.peff")
except PeffParseError as e:
    print(e.line)     # 1-based line number where it failed
    print(e.context)  # the exact offending text
    print(e.hint)     # a short suggestion for how to fix it
except PeffError:
    ...               # any other pefftacular failure
```

`write_peff` validates its input before writing anything and raises
`PeffWriteError` (also a `PeffError`), with a `.hint` and the bad entry's 0-based
`.index` (message prefix `Entry N: `), for a missing header or
an entry with an empty prefix, `db_unique_id` or sequence. File-system failures
(missing directory, no permission) are not wrapped: they raise the usual
`OSError` (`FileNotFoundError`, `PermissionError`).

```python
from pefftacular import PeffWriteError

try:
    write_peff(file_header, entries, "output.peff")
except PeffWriteError as e:  # invalid header/entries
    print(e, e.hint)
except OSError as e:          # could not open/write the file
    print(e)
```

## Spec-violation warnings

Reading is **permissive**: the data is always returned, but anything that
violates a PEFF `MUST` rule (out-of-range positions, missing required fields,
`NumberOfEntries` mismatches, un-coercible custom values, …) is reported through
the `PeffWarning` category. Promote them to errors when you want strict parsing:

```python
import warnings
from pefftacular import PeffWarning, read_peff

warnings.simplefilter("error", PeffWarning)
header, entries = read_peff("suspect.peff")  # now raises on any spec violation
```

## Logging

The library follows the standard logging convention (it attaches a
`NullHandler` and never configures logging itself). Enable a behavioral trace —
useful when scripting or debugging with an AI coding agent:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("pefftacular").setLevel(logging.DEBUG)
```

Milestones (entries read/written) log at `INFO`; file open, header parse, and
entry counts log at `DEBUG`, under the `pefftacular.parser` / `pefftacular.writer`
loggers.

## Development

Contributor and AI-agent guidance lives in [CLAUDE.md](https://github.com/tacular-omics/pefftacular/blob/main/CLAUDE.md). The one
command to run before committing is `just check` (formatting, lint, types, and
tests — the same gate CI enforces); `just fix` auto-applies formatting.

```bash
just install      # install dependencies
just check        # format-check + lint + type-check + test (pre-commit gate)
just fix          # auto-fix lint + formatting
just test         # run tests
just test-file tests/test_errors.py   # run a single test file
just cov          # run tests with coverage
just build        # build the package
just clean        # remove cache files
```

Run `just` with no arguments to list every recipe.

## Citation

If you use pefftacular in research, please cite the archived software release. Machine-readable citation metadata is available in [`CITATION.cff`](https://github.com/tacular-omics/pefftacular/blob/main/CITATION.cff); GitHub's **Cite this repository** menu can render it as APA or BibTeX. DOI: [10.5281/zenodo.22925639](https://doi.org/10.5281/zenodo.22925639).

## License

[MIT](https://github.com/tacular-omics/pefftacular/blob/main/LICENSE)

## Funding

Supported by NIH grants R01AG077046, R01MH132570, R01MH100175, R01HL165168 and U01AG088679.
