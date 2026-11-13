# Tables with pandas or polars

`to_records(source)` returns one plain `dict` per entry, so any data-frame library can
take the result directly. pefftacular does not ship or require pandas or polars;
install whichever you use. (The test suite runs these examples only when the library is
installed.)

```python
import pandas as pd
import polars as pl

from pefftacular import to_records

records = to_records("proteins.peff")   # a path, an open text handle, or entries
df = pd.DataFrame(records)
with_variants = df[df["variant_simple"].notna()]    # entries that have \VariantSimple
human = pl.DataFrame(records).filter(pl.col("ncbi_tax_id") == 9606)
```

`PeffReader.to_records()` and `SequenceEntry.to_record()` give the same dicts. The file
header is not included. Every record has these keys, in this order
(`pefftacular.RECORD_KEYS`); an absent value is `None`:

| key | type | from |
|---|---|---|
| `prefix`, `db_unique_id` | str | the identifier `prefix:db_unique_id` |
| `id`, `db_unique_id_key` | str | `\ID`, `\DbUniqueId` |
| `pname`, `gname`, `tax_name`, `comment` | str | `\PName`, `\GName`, `\TaxName`, `\Comment` |
| `ncbi_tax_id`, `length`, `sv`, `ev`, `pe` | int | `\NcbiTaxId`, `\Length`, `\SV`, `\EV`, `\PE` |
| `decoy` | bool | `\Decoy` |
| `variant_simple`, `variant_complex`, `mod_res_unimod`, `mod_res_psi`, `mod_res`, `processed`, `disulfide_bond`, `proteoform` | str | the annotation's PEFF value text, e.g. `(12\|L)(30\|*)` |
| `custom_values`, `extra` | str | header-declared custom keys / unknown keys as `\Key=value` text |
| `sequence` | str | the residues |

The same fields in [fastatacular](https://github.com/tacular-omics/fastatacular) records have
other names where each package follows its own model. Rename these to put both in one
frame:

| field | fastatacular | pefftacular |
|---|---|---|
| database prefix | `prefix` | `prefix` |
| accession | `accession` | `db_unique_id` |
| entry name | `entry_name` | `id` |
| protein name, gene | `pname`, `gname` | `pname`, `gname` |
| organism name | `os_name` | `tax_name` |
| taxon id, PE, SV | `ncbi_tax_id`, `pe`, `sv` | `ncbi_tax_id`, `pe`, `sv` |
| other keys | `extra`, `KEY=value` pairs | `extra`, `\Key=value` pairs |
| length, residues | `length`, `sequence` | `length`, `sequence` |


## API

::: pefftacular.to_records

::: pefftacular.PeffReader.to_records

::: pefftacular.SequenceEntry.to_record
