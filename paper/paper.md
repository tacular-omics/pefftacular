---
title: 'pefftacular: Reading, writing and validating PSI Extended FASTA (PEFF) files in Python'
tags:
  - Python
  - Proteomics
  - Mass Spectrometry
  - PEFF
  - FASTA
  - Proteoforms
authors:
  - name: Patrick T. Garrett
    orcid: 0000-0002-8434-9693
    affiliation: 1
  - given-names: John R.
    surname: Yates
    suffix: III
    orcid: 0000-0001-5267-1672
    corresponding: true
    affiliation: 1
affiliations:
  - name: The Scripps Research Institute, United States
    index: 1
date: 24 September 2026
bibliography: paper.bib
---

# Summary

The PSI Extended FASTA Format (PEFF) is the HUPO Proteomics Standards Initiative's
standard for protein sequence databases that carry annotations [@binz-2019]. A PEFF file
is a FASTA file with a structured file header and `\Key=value` annotations in each entry
header. The annotations record sequence variants, post-translational modifications,
processed chains, disulfide bonds and proteoforms, so that a search engine can consider
them without writing every variant sequence out as a separate entry. **pefftacular** is a
dependency-free Python library that reads PEFF files into typed, immutable objects, writes
them back, checks them against the specification, and converts entries to and from plain
FASTA and ProForma notation [@leduc-2022].

# Statement of need

PEFF annotations use nested, parenthesised values, backslash escapes, and references
between annotations through `AnnotId` identifiers. A variant list such as
`\VariantSimple=(2|T|7)(3|R|8)...` or a proteoform that points to three other annotations
cannot be handled reliably with the regular expressions that are usually applied to
FASTA headers. The specification also separates rules that a file MUST follow from
recommendations, and it lets each database define its own keys through a header
`CustomKeyDef`, with a regular expression and typed fields.

Groups that build PEFF databases, for example from neXtProt or UniProt exports or from
sample-specific variant calls, need to read those files, change or filter entries, and
write valid PEFF again. Groups that analyse search results need the annotations as data:
which residues can carry which modification, and which variants are possible at a site.
pefftacular gives Python users both directions, with errors that name the line and the
problem, and it targets developers of database-building tools, search pipelines, and
proteoform analyses.

# State of the field

The PSI maintains a list of PEFF implementations [@psi-peff]. Data providers include
neXtProt, which exports PEFF, and the UniProt Proteins API. Several search engines read
PEFF, including Comet [@eng-2020], Protein Prospector and ProteinPilot. Other listed tools
are the Proteomics::PEFF Perl module, the PeptideAtlas online PEFF validator, Proteoformer
and PrecisionProDB. Most of these consume PEFF inside a larger application or a different
language, and none is a general Python library.

In Python, **Pyteomics** [@goloborodko-2013; @levitsky-2019] provides `IndexedPEFF`
[@pyteomics-peff]. It reads a PEFF file with random access and returns each entry's
header as generic key/value pairs with type coercion. It does not write PEFF, and it does
not model the structure of variant, modification or proteoform annotations.
**Biopython** [@cock-2009] reads PEFF only as FASTA, with the whole annotated header kept
as a description string.

pefftacular covers what these do not: a typed model for every annotation key in PEFF 1.0,
typed parsing of header-declared custom keys, a writer that validates its output, and
conversions to ProForma and plain FASTA. Its companion package fastatacular
[@fastatacular] has the same reader and writer API for plain FASTA.

# Software design

**Reading.** `read_peff` loads a whole file and `PeffReader` streams entries one at a
time. Paths may be plain, gzip, bzip2 or xz files. Compression is detected from magic
bytes, so pipes also work. The file header becomes a `FileHeader` with one
`DatabaseHeader` per database section, including `CustomKeyDef` and `OptionalTagDef`
declarations. Each entry becomes a frozen `SequenceEntry`. The standard keys are typed
fields: `PName`, `GName`, `NcbiTaxId`, `Length`, `SV`, `EV` and `PE`. Annotation keys become
tuples of frozen dataclasses:

- `VariantSimple` and `VariantComplex` for substitutions and multi-residue changes;
- `ModResUnimod`, `ModResPsi` and `ModRes` for modification sites, referring to Unimod
  [@creasy-2004], PSI-MOD [@montecchi-palazzi-2008] or other vocabularies;
- `Processed` for chains, signal peptides and other processed forms;
- `DisulfideBond` and `Proteoform`, which refer to other annotations by `AnnotId`.

Values of a key declared by a `CustomKeyDef` are split with its regular expression and
converted to the declared XSD types (string, integer, decimal, boolean, date, time, or an
enumeration). Undeclared keys are kept as text in `extra`, so no information is lost.
Escaped characters are decoded as described in section 3.3.3 of the specification.

**Strict and permissive reading.** Real PEFF files do not always follow every rule. By
default, a violation of a MUST rule, for example a `Length` that disagrees with the
sequence or the deprecated `\Variant=` key, raises a `PeffWarning` and reading continues.
Users who need strict validation turn these warnings into errors with Python's standard
`warnings` filters. Text that cannot be parsed at all raises `PeffParseError`, which
carries the line number, the offending text and a hint for fixing it.

**Writing.** `write_peff` formats and checks every entry before it writes anything, so a
failed write leaves no partial file. By default it also parses each formatted entry back
and compares it with the original. A value that contains PEFF syntax, such as a
`\Key=` token in a protein name, therefore raises `PeffWriteError` instead of producing a
file that reads back differently. The check can be switched off for entries that were
read unchanged from a valid file.

**Conversions.** `SequenceEntry.from_fasta` builds an entry from a UniProt-style FASTA
header, mapping `OS`, `OX`, `GN`, `PE` and `SV` to their PEFF keys, and `to_fasta` goes
back. `to_proforma` writes the sequence with its PSI-MOD or Unimod sites as a ProForma 2.0
string [@leduc-2022], optionally with selected substitutions applied. Both vocabularies
are selected by name, and entries whose sites cannot be written can be skipped instead of
raising an error. `to_records`
returns one flat dictionary per entry for pandas or polars, without depending on either.

**Testing.** The test suite parses every example in the PEFF 1.0 specification and the
example files published by the PSI, including invalid ones, for which it checks the
expected warnings.
Hypothesis property tests [@maciver-2019] check write-read round trips. Continuous
integration runs on Python 3.12 to 3.14 on Linux, and on macOS and Windows, and also
tests the built wheel.

# Research impact statement

pefftacular is part of the tacular-omics packages. Its documentation shows how to turn a
UniProt FASTA file read with fastatacular [@fastatacular] into a PEFF database, and its
ProForma output can be read by Peptacular [@peptacular], for example to compute
the mass of the modified protein.

<!-- TODO(author): add evidence of research use (published analyses, other groups,
lab pipelines, PEFF databases built with it) if any exists. None is recorded in the
repository, and JOSS requires evidence that the software is used for research. -->

# Example usage

The examples use the insulin entry from section 3.4.2 of the PEFF 1.0 specification,
which is included in the repository's test fixtures. The output shown is the real output.

```python
from pefftacular import read_peff

header, [ins] = read_peff("insulin.peff")
print(header.databases[0].prefix, ins.db_unique_id, ins.gname, ins.length)
print(len(ins.variant_simple), ins.variant_simple[0])
print(ins.mod_res_psi[0]); print(ins.processed[1]); print(ins.disulfide_bond[0])
print(ins.proteoform[10].ranges, ins.proteoform[10].annot_id_refs)
```

```text
nxp NX_P01308-1 INS 110
70 VariantSimple(position=2, new_amino_acid='T', tag=None, annot_id=7)
ModResPsi(positions=(53,), accession='MOD:00087', name='N6-myristoyl-L-lysine', tag=None, annot_id=0)
Processed(start_pos=25, end_pos=54, accession='PEFF:0001020', name='mature protein', tag=None, annot_id=78)
DisulfideBond(annot_id_refs=(1, 2), description='between chains', annot_id=81)
(SequenceRange(start=90, end=110), SequenceRange(start=25, end=54)) (81, 82, 83)
```

```python
import warnings
from pefftacular import PeffWarning, SequenceEntry, read_peff, write_peff

header, [ins] = read_peff("insulin.peff")
print(ins.to_proforma()[:44])
print(ins.to_proforma(variants=ins.variant_simple[:1])[:12])

sp = SequenceEntry.from_fasta(
    "sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1",
    ins.sequence,
)
print(sp.prefix, sp.db_unique_id, sp.id, sp.ncbi_tax_id)

write_peff(header, [ins], "copy.peff")          # each entry is read back and checked
with warnings.catch_warnings():
    warnings.simplefilter("error", PeffWarning)  # strict mode
    assert read_peff("copy.peff")[1] == [ins]
```

```text
MALWMRLLPLLALLALWGPDPAAAFVNQHLC[MOD:00798]GS
MTLWMRLLPLLA
sp P01308 INS_HUMAN 9606
```

# Availability

pefftacular is released under the MIT license. It is distributed on PyPI
(<https://pypi.org/project/pefftacular/>), developed on GitHub
(<https://github.com/tacular-omics/pefftacular>), documented at
<https://tacular-omics.github.io/pefftacular/>, and archived on Zenodo [@pefftacular]. It
requires Python 3.12 or later.

# AI usage disclosure

During the preparation of this work the authors used Anthropic Claude large language models via the Claude Code interface for software-development assistance, including code, tests, and documentation, and for manuscript drafting and editing. The authors reviewed and edited all content and take full responsibility for the software and the publication.

# Acknowledgements

We thank Claire Delahunty, Ph.D., for a careful reading of the manuscript. This work was supported by the U.S. National Institutes of Health (grants R01 HL165168, R01 AG077046, R01 MH100175, and R01 AG075862 to J.R.Y.) and by the Skaggs Graduate School of Chemical and Biological Sciences at The Scripps Research Institute (P.T.G.).

# References
