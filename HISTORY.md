# History

## 0.4.0 (2026-07-10)

### Breaking

* `DisulfideBond.positions` renamed to `DisulfideBond.annot_id_refs` to reflect that a `\DisulfideBond` references prior `ModResPsi` entries by annotation ID, not by residue position (spec section 3.4.2).

### Spec compliance

* Entry-item components now honor the spec's backslash escaping (section 3.3.3): `\|`, `\(`, `\)`, and `\\` parse as literals, and `write_peff()` emits those escapes so values containing pipes or unpaired parens round-trip. Balanced parens (e.g. `N-linked (GlcNAc...)`) are left unescaped.
* Database-level `# GeneralComment=` lines are now preserved on `DatabaseHeader.general_comments` and re-emitted, instead of being silently dropped.
* `DbName` is now included in the missing-mandatory-key warning for database headers.
* `ProteoformDb` and `HasAnnotationIdentifiers` header flags are now recognized case-insensitively (tolerating the spec's own `ProteoformDB` / singular-form variants).

### Errors, warnings, and logging

* Added a `PeffError` base class (subclass of `ValueError`); `PeffParseError` and `PeffWriteError` now derive from it, so any library failure can be caught with a single `except PeffError`.
* Parse/write errors now carry an actionable `hint`, and `PeffParseError` attaches its `context` and `hint` as exception notes (PEP 678) so they appear in tracebacks.
* Non-fatal spec violations are now emitted through a dedicated `PeffWarning` category (subclass of `UserWarning`), so they can be filtered or escalated on their own: annotation `MUST`-rule checks (positions outside `1..len(sequence)`, empty `VariantSimple` `newAminoAcid`, missing required accession/name on `ModRes*`) plus the existing version/header/count/custom-value warnings.
* Added `logging` throughout the reader and writer under the `pefftacular.parser` / `pefftacular.writer` loggers (with a package-level `NullHandler`), for a behavioral trace at `DEBUG`/`INFO`.

### Tooling & docs

* Added `AGENTS.md` and `CLAUDE.md` guidance for AI coding agents.
* Reworked the `justfile`: `just check` is a single read-only pre-commit gate (format-check + lint + types + tests over `src` and `tests`), and `just fix` auto-applies fixes.
* Documented error handling, spec-violation warnings, and logging in the README.

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
