# Annotations

PEFF's value over FASTA is the per-entry annotation: sequence variants, modifications,
processing events, disulfide bonds and proteoforms. pefftacular parses each annotation key into a
tuple of typed objects on the `SequenceEntry`.

| PEFF key | Field on `SequenceEntry` | Model | Describes |
|---|---|---|---|
| `\VariantSimple` | `variant_simple` | `VariantSimple` | a single-residue substitution |
| `\VariantComplex` | `variant_complex` | `VariantComplex` | an insertion, deletion or multi-residue replacement |
| `\ModResUnimod` | `mod_res_unimod` | `ModResUnimod` | a modification with a UNIMOD accession |
| `\ModResPsi` | `mod_res_psi` | `ModResPsi` | a modification with a PSI-MOD accession |
| `\ModRes` | `mod_res` | `ModRes` | a modification with any other (or no) accession |
| `\Processed` | `processed` | `Processed` | a processed region: signal peptide, chain, propeptide |
| `\DisulfideBond` | `disulfide_bond` | `DisulfideBond` | a bond between two annotated half-cystines |
| `\Proteoform` | `proteoform` | `Proteoform` | a defined combination of ranges and annotations |
| declared by `CustomKeyDef` | `custom_values` | `CustomKeyValue` | your own structured keys |

## The example entry

Every example on this page reads the same insulin entry, adapted from the official PEFF
examples. Run this first:

```python
import io

from pefftacular import read_peff

PEFF = r"""# PEFF 1.0
# //
# DbName=Insulin example
# Prefix=nxp
# DbVersion=1.0
# DbSource=www.nextprot.org
# NumberOfEntries=1
# SequenceType=AA
# HasAnnotationIdentifiers=true
# //
>nxp:NX_P01308-1 \PName=Insulin \GName=INS \NcbiTaxId=9606 \Length=110 \ModResPsi=(1:31|MOD:00798|half cystine)(2:96|MOD:00798|half cystine) \ModResUnimod=(3:53|UNIMOD:45|Myristoyl) \VariantSimple=(4:6|C)(5:12|V|dbSNP) \VariantComplex=(6:1|3|M) \Processed=(7:1|24|PEFF:0001021|signal peptide)(8:25|54|PEFF:0001020|mature protein) \DisulfideBond=(9:1,2|between chains) \Proteoform=(10:NX_P01308-1-pf1|1-110||preproinsulin)(11:NX_P01308-1-pf2|90-110,25-54|9|insulin A and B chains joined)
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
"""

header, (insulin,) = read_peff(io.StringIO(PEFF))
print(insulin.gname, len(insulin.sequence))
# INS 110
```

## Positions and annotation IDs

Two rules apply to every annotation type:

- **Positions are 1-based and inclusive**, counted on `entry.sequence`. Position 1 is the first
  residue. A position the source marks as unknown (`?`) is kept as the string `"?"`, so position
  fields are typed `int | str`.
- **Annotation IDs** are optional integers written as `id:` in front of the first field
  (`(4:6|C)` is annotation 4, at position 6). They appear when the database header sets
  `HasAnnotationIdentifiers=true`, and are how `DisulfideBond` and `Proteoform` refer to other
  annotations. Without them, `annot_id` is `None`.

```python
print(header.databases[0].has_annotation_identifiers)
# True
print([v.annot_id for v in insulin.variant_simple])
# [4, 5]
```

Each annotation also has an optional `tag`, a free-text source label such as `dbSNP`. If the
database header declares tags with `# OptionalTagDef=tag:description`, they are available as
`DatabaseHeader.optional_tag_defs`.

## Variants

### `VariantSimple`: one residue changes

`(position|newAminoAcid|tag)`. `new_amino_acid` is a single residue, or `*` for a stop codon.

```python
for v in insulin.variant_simple:
    ref = insulin.sequence[v.position - 1]
    print(f"{ref}{v.position}{v.new_amino_acid}", v.tag)
# R6C None
# A12V dbSNP
```

`apply_variants` builds the variant sequence:

```python
form = insulin.apply_variants([insulin.variant_simple[1]])
print(form.sequence[:15], len(form.sequence))
# MALWMRLLPLLVLLA 110
```

A `VariantSimple` never changes the length (spec §3.3.8). A `*` is a nonsense mutation: the
sequence ends before that position and `form.truncated_at` holds it.

### `VariantComplex`: a stretch changes

`(startPosition|endPosition|newSequence|tag)`. Residues `start_pos` to `end_pos` inclusive are
replaced by `new_sequence`. An empty `new_sequence` is a deletion; a `new_sequence` longer than the
range is an insertion.

```python
vc = insulin.variant_complex[0]
print(vc.start_pos, vc.end_pos, repr(vc.new_sequence))
# 1 3 'M'

form = insulin.apply_variants([vc])
print(form.sequence[:10], len(form.sequence))
# MWMRLLPLLA 108
```

### Applying variants

`apply_variants` returns a `VariantForm`. Its `entry` is an ordinary `SequenceEntry`, so
`to_proforma`, `to_record` and `write_peff` work on it unchanged, and `position_map` says where
every original residue went: 1-based on both sides, `None` where the variant removed it.

```python
print(form.position_map[1], form.position_map[3], form.position_map[4])
# None None 2
```

Every other annotation is carried onto the new coordinates, and each one that cannot be is
dropped with a `PeffWarning` rather than left pointing at the wrong residue:

```python
print([(p.name, p.start_pos, p.end_pos) for p in form.entry.processed])
# [('mature protein', 23, 52)]
print([(m.positions, m.name) for m in form.entry.mod_res_psi])
# [((29,), 'half cystine'), ((94,), 'half cystine')]
print(form.entry.length)
# 108
```

The signal peptide (1-24) and the `preproinsulin` proteoform (1-110) are gone: both start at
position 1, which this variant replaced, and a range that lost a bound is dropped rather than
quietly clipped. A modification on a residue a `VariantSimple` substituted is dropped too, because
the spec says such a proteoform needs its own entry (§3.3.10).

### Enumerating the forms

A variant is optional unless the database says otherwise ("MAY be substituted", §3.3.8), so
`variant_forms` yields the entry as written and then one form per variant:

```python
forms = list(insulin.variant_forms())
print(len(forms), forms[0].is_reference)
# 4 True
print([len(f.applied) for f in forms])
# [0, 1, 1, 1]
```

Combinations grow fast, so they are opt-in and bounded. `max_combination=k` applies up to `k`
variants at once, `limit` (10,000 by default) caps the count and is checked before a single
sequence is built, and overlapping variants are skipped with a warning, since the spec does not
define how they combine. For the forms the file itself means, pass `proteoforms=True`, which
applies the variants each `\Proteoform` references by annotation id (§3.4.2):

```python
print(len(list(insulin.variant_forms(max_combination=3))))
# 8
print(len(list(insulin.variant_forms(proteoforms=True))))
# 2
```

The difference matters on real files. The spec's own insulin example
(`tests/fixtures/PEFF_AnnotID_Insulin_Valid.peff`) carries 70 `VariantSimple` annotations: 71
forms by default, 2,486 for every pair, 2\*\*70 for the full power set, and 11 proteoforms it
declares itself. For a database whose header has `ProteoformDb=true`, where every annotation MUST
be applied (§3.4.1), pass `mandatory=True` to get that one sequence.

## Modifications

The three modification keys share one shape, `(positions|accession|name|tag)`, and differ only in
the controlled vocabulary of the accession:

- `ModResUnimod`: a UNIMOD accession, such as `UNIMOD:21` (Phospho).
- `ModResPsi`: a PSI-MOD accession, such as `MOD:00798` (half cystine).
- `ModRes`: any other accession, or none. Use it for modifications that are not in either
  vocabulary.

`positions` is always a tuple, because one annotation can mark the same modification at several
residues: `(12,15,20|UNIMOD:21|Phospho)`.

```python
for mod in insulin.mod_res_unimod:
    print(mod.positions, mod.accession, mod.name)
# (53,) UNIMOD:45 Myristoyl

for mod in insulin.mod_res_psi:
    print(mod.annot_id, mod.positions, mod.accession, mod.name)
# 1 (31,) MOD:00798 half cystine
# 2 (96,) MOD:00798 half cystine
```

To collect every modified position, whatever its vocabulary:

```python
sites = sorted(
    (p, m.name)
    for m in (*insulin.mod_res_unimod, *insulin.mod_res_psi, *insulin.mod_res)
    for p in m.positions
)
print(sites)
# [(31, 'half cystine'), (53, 'Myristoyl'), (96, 'half cystine')]
```

## Processed regions

`Processed` marks a region produced by post-translational processing:
`(startPosition|endPosition|accession|name|tag)`. The accession comes from the PEFF controlled
vocabulary, for example `PEFF:0001021` (signal peptide) and `PEFF:0001020` (mature protein).

```python
for region in insulin.processed:
    piece = insulin.sequence[region.start_pos - 1 : region.end_pos]
    print(region.name, region.start_pos, region.end_pos, piece[:8])
# signal peptide 1 24 MALWMRLL
# mature protein 25 54 FVNQHLCG
```

## Disulfide bonds

A `\DisulfideBond` does not give residue positions. It references two earlier `ModResPsi`
half-cystine annotations **by annotation ID**, so `annot_id_refs` holds IDs, not positions.
Resolve them through the modifications:

```python
by_id = {m.annot_id: m for m in insulin.mod_res_psi}

for bond in insulin.disulfide_bond:
    a, b = (by_id[ref].positions[0] for ref in bond.annot_id_refs)
    print(f"Cys{a}-Cys{b}", bond.description)
# Cys31-Cys96 between chains
```

## Proteoforms

A `Proteoform` names one specific molecular form of the protein:
`(proteoformId|ranges|annotIdRefs|name)`.

- `ranges` is a tuple of `SequenceRange(start, end)`. Several ranges mean several chains, such as
  insulin's A and B chains.
- `annot_id_refs` lists the annotation IDs (modifications, variants, bonds) present in this form.

```python
for pf in insulin.proteoform:
    spans = [(r.start, r.end) for r in pf.ranges]
    print(pf.proteoform_id, spans, pf.annot_id_refs, pf.name)
# NX_P01308-1-pf1 [(1, 110)] () preproinsulin
# NX_P01308-1-pf2 [(90, 110), (25, 54)] (9,) insulin A and B chains joined
```

Databases built around proteoforms set `ProteoformDb=true` in the header, available as
`DatabaseHeader.proteoform_db`.

## Custom keys

A database header can declare its own entry keys with `# CustomKeyDef=`. The definition names
the fields, their XSD types, and optionally a regular expression that splits the value.
pefftacular uses it to parse matching entry values into typed fields on `entry.custom_values`:

```python
import io

from pefftacular import read_peff

CUSTOM = r"""# PEFF 1.0
# //
# DbName=Custom demo
# Prefix=cu
# DbVersion=1
# DbSource=local
# NumberOfEntries=1
# SequenceType=AA
# CustomKeyDef=(KeyName=SecondaryStructure|Description="Secondary structure element"|RegExp="([0-9]+)\|([0-9]+)\|(.+)"|FieldNames=StartPosition,EndPosition,Description|FieldTypes=integer,integer,string)
# //
>cu:P00001 \SecondaryStructure=(10|20|Helix)(25|30|Strand)
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK
"""

header, (entry,) = read_peff(io.StringIO(CUSTOM))

key_def = header.databases[0].custom_key_defs[0]
print(key_def.key_name, key_def.field_names, key_def.field_types)
# SecondaryStructure ('StartPosition', 'EndPosition', 'Description') ('integer', 'integer', 'string')

for value in entry.custom_values["SecondaryStructure"]:
    print(value.fields, value.raw)
# {'StartPosition': 10, 'EndPosition': 20, 'Description': 'Helix'} 10|20|Helix
# {'StartPosition': 25, 'EndPosition': 30, 'Description': 'Strand'} 25|30|Strand
```

Supported `FieldTypes` are `string`, `integer`, `decimal`, `boolean`, `date`, `time` and
`enumeration(a|b|c)`. A value that does not fit its type stays a string and raises a
`PeffWarning`: a `boolean` field accepts `true`, `1`, `yes`, `false`, `0` and `no`,
case-insensitively, and any other value warns and stays a string rather than being read as
`False`. With no `RegExp`, the value is split on `|` and paired with `FieldNames` in order.

`raw` keeps the original text. The writer uses it only while it still matches `fields`, so a
read-then-write round trip reproduces the value exactly, and an edit made with
`dataclasses.replace(value, fields=...)` is written. To write a value you built yourself, leave
`raw` empty and fill `fields`.

Keys that no `CustomKeyDef` declares are not dropped: they land in `entry.extra` as raw strings.

```python
print(entry.extra)
# {}
```

## ProForma strings

`entry.to_proforma()` writes the sequence with its modification sites as a
[ProForma 2.0](https://github.com/HUPO-PSI/ProForma) string. Pass it to a ProForma parser such as
[peptacular](https://github.com/tacular-omics/peptacular) to get masses or fragments.

```python
print(insulin.to_proforma()[25:50])
# VNQHLC[MOD:00798]GSHLVEAL
print(insulin.to_proforma(mods="unimod")[44:64])
# ERGFFYTPK[UNIMOD:45]
```

- `mods="psimod"` (the default) writes `\ModResPsi` sites and `mods="unimod"` writes
  `\ModResUnimod` sites. A `\ModRes` site is included when its accession is from the same
  vocabulary (`MOD:` or `UNIMOD:`), and left out otherwise.
- Every listed site is modified at once. To render a subset, build a copy with
  `dataclasses.replace(entry, mod_res_psi=...)` first.
- A `?` position becomes a ProForma unknown-position modification: `[MOD:00046]^2?SEQ...`.
- A site listed in both `\ModResPsi` (or `\ModResUnimod`) and `\ModRes` is written once.
- A modification with no accession is written by name: `[M:name]` or `[U:name]`.
- PEFF does not say whether a modification on residue 1 is on the N-terminus or on the side
  chain, so it is always written on the residue.

`variants=` applies `VariantSimple` substitutions first. A modification on a substituted
residue is dropped (the spec says a modified variant needs its own entry), and a `*` variant
truncates the sequence before its position:

```python
print(insulin.to_proforma(variants=insulin.variant_simple[:1])[:12])
# MALWMCLLPLLA
```

`?` sites are dropped after a `*` truncation, since they may lie in the removed part.

Positions outside the sequence, or two different substitutions at one position, raise
`PeffError`; the message starts with the entry's `prefix:db_unique_id`. `VariantComplex`
is not applied. Real files have such entries (12 of the 20,431 in the neXtProt human
PEFF list sites past the sequence end). To convert a whole file, skip them with
`errors="skip"`, which returns `None` for an entry that cannot be written:

```python
forms = [p for e in entries if (p := e.to_proforma(errors="skip")) is not None]
```

## Writing annotations

Build the annotation objects and pass them to `SequenceEntry` as tuples:

```python
import io

from pefftacular import (
    DatabaseHeader,
    FileHeader,
    ModResPsi,
    Processed,
    SequenceEntry,
    VariantSimple,
    write_peff,
)

header = FileHeader(
    peff_version="1.0",
    databases=(DatabaseHeader(prefix="my", db_name="Demo", db_version="1", db_sources=("local",),
                              number_of_entries=1, sequence_type="AA"),),
)

entry = SequenceEntry(
    prefix="my",
    db_unique_id="PROT001",
    sequence="MKTAYIAKQRQISFVKSHFSRQ",
    variant_simple=(VariantSimple(position=5, new_amino_acid="C"),),
    mod_res_psi=(ModResPsi(positions=(13, 16), accession="MOD:00046", name="O-phospho-L-serine"),),
    processed=(Processed(start_pos=1, end_pos=4, accession="PEFF:0001021", name="signal peptide"),),
)

out = io.StringIO()
write_peff(header, [entry], out)
print(out.getvalue().splitlines()[-2])
# >my:PROT001 \VariantSimple=(5|C) \ModResPsi=(13,16|MOD:00046|O-phospho-L-serine) \Processed=(1|4|PEFF:0001021|signal peptide)
```

Set `annot_id` on each annotation if you want annotation IDs in the output, and set
`HasAnnotationIdentifiers` by passing `has_annotation_identifiers=True` to the `DatabaseHeader`.

The writer escapes `\` and `|` inside annotation fields, and parentheses when they are
unbalanced, so names such as `N-linked (GlcNAc...)` or `a|b` survive a round trip.

The same escaping applies to the plain-text fields `pname`, `gname`, `tax_name` and
`comment`, and to custom-key fields built without a `RegExp`; the reader reverses it.

!!! warning "`extra` values are written as-is"
    Values in `extra` are raw strings and are written without escaping. Keep them free of
    the sequence ` \` (space, backslash), which starts a new key.
