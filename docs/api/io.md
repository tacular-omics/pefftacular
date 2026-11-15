# Reading and writing

::: pefftacular.read_peff

::: pefftacular.PeffReader

::: pefftacular.write_peff

## Records

`to_records` turns entries into flat dicts for a data frame (see the README section on
tables). `from_records` is its inverse and takes the file header, whose `CustomKeyDef`
declarations are what keep custom key values from being read back as unknown keys.

::: pefftacular.from_records

## Compression

`Compression` types the `compression` argument of `read_peff`, `PeffReader`, `write_peff`
and `to_records`.

::: pefftacular.Compression
