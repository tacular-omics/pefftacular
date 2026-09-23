# From FASTA and UniProt

pefftacular has no built-in converter: it reads and writes PEFF, and does not fetch or parse
other formats. Converting is short, though, because a PEFF entry is just a `SequenceEntry` you
build yourself. This page shows two ways to do it.

## From a UniProt FASTA file with fastatacular

[fastatacular](https://github.com/tacular-omics/fastatacular) is the plain-FASTA sibling of
pefftacular. It parses UniProt-style headers (`>sp|P69905|HBA_HUMAN ... OS=... OX=... GN=...`)
into named fields, which map almost one-to-one onto PEFF keys. Install it alongside
pefftacular:

```bash
pip install pefftacular fastatacular
```

Take this UniProt FASTA file, saved as `uniprot.fasta`:

```text
>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH
```

Map each FASTA record to a `SequenceEntry` and write the PEFF file:

```python
from fastatacular import read_fasta
from pefftacular import DatabaseHeader, FileHeader, SequenceEntry, write_peff

records = read_fasta("uniprot.fasta")

entries = [
    SequenceEntry(
        prefix="sp",
        db_unique_id=r.accession,
        sequence=r.sequence,
        pname=r.pname,
        gname=r.gname,
        ncbi_tax_id=r.ncbi_tax_id,
        tax_name=r.os_name,
        length=len(r.sequence),
        sv=r.sv,
        pe=r.pe,
    )
    for r in records
]

header = FileHeader(
    peff_version="1.0",
    databases=(
        DatabaseHeader(
            prefix="sp",
            db_name="UniProtKB/Swiss-Prot",
            db_version="2026_03",
            db_sources=("https://www.uniprot.org",),
            number_of_entries=len(entries),
            sequence_type="AA",
        ),
    ),
)

write_peff(header, entries, "uniprot.peff")
print(open("uniprot.peff").read().splitlines()[9])
# >sp:P69905 \Length=142 \PName=Hemoglobin subunit alpha \GName=HBA1 \NcbiTaxId=9606 \TaxName=Homo sapiens \SV=2 \PE=1
```

A few details worth knowing:

- UniProt FASTA files mix Swiss-Prot (`sp`) and TrEMBL (`tr`) records. `r.prefix` holds the
  source, so you can use `prefix=r.prefix` and write one `DatabaseHeader` per prefix.
- `NumberOfEntries` must match the number of entries per prefix, or readers will warn.
- UniProt FASTA carries no variants or modifications. Add them from another source (UniProt's
  feature table, dbSNP, your own search results) by building `VariantSimple`, `ModResUnimod` and
  other annotation objects, as shown in the [annotations guide](annotations.md#writing-annotations).

## From any FASTA without extra dependencies

If you do not want another package, a FASTA reader is a few lines. This version keeps the first
word of each header as the ID and the rest as the protein name:

```python
from pefftacular import DatabaseHeader, FileHeader, SequenceEntry, write_peff


def read_simple_fasta(path):
    name, desc, seq = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(seq)
                name, _, desc = line[1:].partition(" ")
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, desc, "".join(seq)


entries = [
    SequenceEntry(
        prefix="lab",
        db_unique_id=name.replace("|", "_"),
        sequence=seq,
        pname=desc or None,
        length=len(seq),
    )
    for name, desc, seq in read_simple_fasta("uniprot.fasta")
]

header = FileHeader(
    peff_version="1.0",
    databases=(
        DatabaseHeader(
            prefix="lab",
            db_name="Lab database",
            db_version="1",
            db_sources=("uniprot.fasta",),
            number_of_entries=len(entries),
            sequence_type="AA",
        ),
    ),
)

write_peff(header, entries, "lab.peff")
print(entries[0].db_unique_id, entries[0].length)
# sp_P69905_HBA_HUMAN 142
```

The part of a PEFF description line after the `prefix:` and before the first space is the
unique ID, so it must not contain spaces. The example replaces `|` too, to keep IDs readable.

## Files that are already PEFF

Some sources publish PEFF directly. The official PEFF example files, for instance, include a
UniProt export that carries `\VariantSimple` annotations and uses UniProt's `\OX` key for the
taxonomy ID (pefftacular maps `\OX` to `ncbi_tax_id`). If you already have a PEFF file, skip
conversion and read it with `read_peff`.
