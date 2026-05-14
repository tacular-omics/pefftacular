# History

## 0.3.0 (2026-05-14)

* Header-declared custom keys (`# CustomKeyDef=`) now drive typed parsing of entry values: registered keys are parsed via their `RegExp` (with pipe-split as fallback) and coerced per `FieldTypes` (`integer`, `decimal`, `boolean`, `date`, `time`, `string`, `enumeration`).
* Parsed custom-key values are exposed on `SequenceEntry.custom_values`, with the raw string preserved alongside for lossless round-trip writing.
* Undeclared custom keys continue to land in `SequenceEntry.extra` unchanged.
* `CustomKeyDef` now captures the `ConceptCURIE` field.
* Multiple `# CustomKeyDef=` lines per database are preserved instead of silently overwriting each other.
* Field tokenization is now quote-aware so a `RegExp` value containing escaped pipes survives splitting.

## 0.2.0 (2026-03-18)

* Added `DisulfideBond`, `Proteoform`, `SequenceRange`, and `OptionalTagDef` models.
* `\DisulfideBond=` and `\Proteoform=` entry keys now parse into typed fields instead of `extra`.
* `\OX=` (UniProt NcbiTaxId alias) now maps to `ncbi_tax_id`.
* `\ID=`, `\DbUniqueId=`, and `\Comment=` entry keys now parse into named fields.
* `\Variant=` (deprecated since 2015) now emits a `DeprecationWarning`.
* `DatabaseHeader` gains named fields: `db_description`, `db_date`, `db_sources`, `decoy`, `conversion`, `has_annotation_identifiers`, `proteoform_db`, `optional_tag_defs`.
* `db_source: str | None` replaced by `db_sources: tuple[str, ...]` to support multiple `DbSource` lines.
* Multi-value header keys (`DbSource`, `OptionalTagDef`) no longer silently overwrite on repeated lines.
* `annot_id: int | None` field added to all annotation models (`VariantSimple`, `VariantComplex`, `ModResUnimod`, `ModResPsi`, `ModRes`, `Processed`) for `HasAnnotationIdentifiers=true` databases.
* Integer parsing errors (`NcbiTaxId`, `Length`, `SV`, `EV`, `PE`, `NumberOfEntries`) now raise `PeffParseError`.
* `Length` mismatch between tag and actual sequence now emits a `UserWarning`.
* `write_peff()` now raises `PeffWriteError` for `None` header or entries with empty `prefix`, `db_unique_id`, or `sequence`.
* Added GitHub Actions CI workflow (Python 3.12 and 3.13).
* Added official PEFF example files as test fixtures.

## 0.1.0 (2026-03-18)

* First release on PyPI.
