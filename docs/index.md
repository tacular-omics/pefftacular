# pefftacular

[![CI](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml/badge.svg)](https://github.com/tacular-omics/pefftacular/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pefftacular)](https://pypi.org/project/pefftacular/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22925639.svg)](https://doi.org/10.5281/zenodo.22925639)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-g.svg)](https://opensource.org/licenses/MIT)

pefftacular reads and writes [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format)
files in pure Python, with no dependencies.

## What is PEFF?

PEFF is the HUPO Proteomics Standards Initiative's extension of FASTA. A PEFF file is still a
FASTA file: every entry is a `>` description line followed by the sequence. PEFF adds two things:

- a **file header** of `# Key=value` lines that describes each source database (name, version,
  prefix, entry count, custom keys), and
- **structured `\Key=value` pairs** on each description line, so an entry can carry gene and
  protein names, taxonomy, sequence variants, post-translational modifications, signal peptides,
  disulfide bonds and whole proteoforms.

```text
# PEFF 1.0
# //
# DbName=SwissProt
# Prefix=sp
# DbVersion=2026_03
# DbSource=https://www.uniprot.org
# NumberOfEntries=1
# SequenceType=AA
# //
>sp:P12345 \PName=Example protein \GName=EXMP \VariantSimple=(5|C|dbSNP) \ModResUnimod=(13|UNIMOD:21|Phospho)
MKTAYIAKQRQISFVKSHFSRQ
```

Search engines that support PEFF, such as Comet, can read these annotations and search variant and
modified peptides directly from the database.

## Why pefftacular?

- **Typed data, not strings.** Headers, entries and every annotation type are frozen dataclasses.
  You get `entry.mod_res_unimod[0].positions`, not a regex over the header.
- **Streaming or eager.** `read_peff` loads a whole file; `PeffReader` yields one entry at a time.
- **Permissive reading, strict when you want it.** Spec violations are reported as
  `PeffWarning` and the data is still returned. Structural failures raise `PeffParseError` with
  the line number, the offending text and a hint.
- **Round-trip writing** with the backslash escaping PEFF requires.
- **Zero dependencies**, Python 3.12+.

## Install

```bash
pip install pefftacular
# or, in a uv project
uv add pefftacular
```

## 30-second example

Build an entry with a variant and a phosphorylation site, write it, and read it back:

```python
from pefftacular import (
    DatabaseHeader,
    FileHeader,
    ModResUnimod,
    SequenceEntry,
    VariantSimple,
    read_peff,
    write_peff,
)

header = FileHeader(
    peff_version="1.0",
    databases=(
        DatabaseHeader(
            prefix="sp",
            db_name="SwissProt",
            db_version="2026_03",
            db_sources=("https://www.uniprot.org",),
            number_of_entries=1,
            sequence_type="AA",
        ),
    ),
)

entry = SequenceEntry(
    prefix="sp",
    db_unique_id="P12345",
    sequence="MKTAYIAKQRQISFVKSHFSRQ",
    pname="Example protein",
    gname="EXMP",
    variant_simple=(VariantSimple(position=5, new_amino_acid="C", tag="dbSNP"),),
    mod_res_unimod=(ModResUnimod(positions=(13,), accession="UNIMOD:21", name="Phospho"),),
)

write_peff(header, [entry], "example.peff")

header, entries = read_peff("example.peff")
for e in entries:
    print(e.db_unique_id, e.gname, len(e.sequence))
    for mod in e.mod_res_unimod:
        print(mod.name, "at", mod.positions)
# P12345 EXMP 22
# Phospho at (13,)

print(entries[0] == entry)
# True
```

## Where next

- [Getting started](getting-started.md): reading, writing, the header and entry model, errors.
- [Annotations](annotations.md): variants, modifications, processing, disulfide bonds,
  proteoforms and custom keys.
- [From FASTA and UniProt](fasta.md): turn a plain FASTA file into PEFF.
- [API reference](api/io.md): every public class and function.
- Using an AI coding assistant? Point it at
  [`llms.txt`](https://github.com/tacular-omics/pefftacular/blob/main/llms.txt) for a compact API guide.

## Related packages

- [fastatacular](https://github.com/tacular-omics/fastatacular) reads and writes plain FASTA with
  the same API shape.
- [peptacular](https://github.com/tacular-omics/peptacular) handles peptide sequences, ProForma
  notation, digestion and mass calculation.
- [unimodpy](https://github.com/tacular-omics/unimodpy) and
  [psimodpy](https://github.com/tacular-omics/psimodpy) look up the UNIMOD and PSI-MOD
  accessions used in `ModResUnimod` and `ModResPsi`.
