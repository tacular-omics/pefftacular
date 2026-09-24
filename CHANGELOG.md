# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Fixed

- `write_peff()` no longer ignores edits to a parsed custom-key value: `dataclasses.replace(value, fields=...)` used to write the stale original text from `CustomKeyValue.raw`. `raw` is now used verbatim only if re-parsing it with the key's `CustomKeyDef` gives the value's current `fields`, otherwise the item is rebuilt from `fields`. Unedited values still round-trip byte-exact.
- The last database block in the file header is no longer dropped when the header runs straight into the first `>` entry without a closing `# //` (now a `PeffWarning`), or when a `# //` separator has trailing whitespace. Previously its `Prefix` and `CustomKeyDef`s were lost and its custom keys ended up in `extra`.
- A blank line inside the file header no longer ends the header and discards the database blocks after it. It is skipped with a `PeffWarning` (spec section 3.3.1: every header line starts with `# `). Blank lines between the header and the first entry are still ignored silently.
- `write_peff()` now escapes a `Proteoform` id, so an id containing `|` or an unbalanced paren no longer writes an unparseable line.
- `write_peff()` no longer writes a trailing `|` for a `DisulfideBond` with `description=""` (it read back as `None`, so the text changed on a second write).
- A quoted `Description` or `RegExp` in a `CustomKeyDef` may now contain an unbalanced paren; it used to raise `PeffParseError` on read.
- `write_peff()` raises `PeffWriteError` instead of writing a corrupt file when an edited field of a `RegExp`-controlled custom key can no longer be read back through the RegExp, when any header or entry value contains a line break, when a prefix contains `:` or whitespace, when a `db_unique_id` contains whitespace, or when a sequence contains whitespace or `>`. The check runs before anything is written.
- The reader now warns (`PeffWarning`) on the spec's illegal examples it used to accept silently: a `VariantSimple` new residue that is not one letter or `*` (section 3.3.8), a `VariantComplex` new sequence with non-residue characters, or a single-residue substitution that should be a `VariantSimple` (3.3.9), and a `Processed` item without an accession or name (3.3.13).

### Tests

- `tests/test_spec_examples.py`: every example in the PEFF 1.0 specification (header, custom keys, entry rules, TYRO3, insulin, legal and illegal annotation examples) parses to the expected structure and writes back; fixtures in `tests/fixtures/spec/`.
- `tests/test_properties.py`: Hypothesis properties for write -> read -> write stability, edited custom-key fields, and malformed input raising only `PeffError`. `HYPOTHESIS_PROFILE=thorough` runs 5000 examples per property.

## [0.4.4] (2026-09-23)

### Fixed

- The description-line key scan now honours spec escapes (section 3.3.3): an escaped unbalanced paren (e.g. a `ModRes` named `odd ( name`, written as `odd \( name`) no longer swallows the following `\Key=`, so such entries round-trip.
- `write_peff()` now escapes the free-text fields `PName`, `GName`, `TaxName` and `Comment` (`\`, `|`, unbalanced parens) and the reader unescapes them, so values containing those characters or ` \Key=` round-trip.
- The deprecated `\Variant=` key now warns with `PeffWarning` instead of `DeprecationWarning`, like every other spec issue.
- `write_peff()` now escapes custom-key fields built from `fields` (keys without a `RegExp`), so a `|` or unbalanced paren in a field no longer reads back as an extra field or a parse error.
- `PeffParseError` raised while parsing an annotation value now carries the entry's line number (`.line`) instead of `None`.
- The `Proteoform` parse-error hint now shows `annotIdRefs` as optional, matching the parser (2 fields accepted).
- README: `ModResUnimod` / `ModResPsi` examples use `positions`, not `position`.
- `CITATION.cff`: `date-released` for 0.4.3 is 2026-09-23.

## [0.4.3] (2026-09-23)

### Changed

- First release archived on Zenodo; no code changes.

## [0.4.2] (2026-09-23)

### Fixed

- Zenodo archiving: removed the grant ids and hard-coded version from `.zenodo.json`, which made Zenodo reject the previous release. Funding is now credited in the README.

## [0.4.1] (2026-09-23)

* Publish from GitHub Actions with PyPI trusted publishing (`publish.yml`);
  release metadata is checked against the tag.
* Keep `__version__` and `CITATION.cff` in sync with
  `scripts/release_version.py` (`just set-version X.Y.Z`).
* CI tests Python 3.12-3.14 on Linux plus macOS and Windows, the lowest
  direct dependency versions, and the built wheel.
* Standardize citation and package metadata.

## [0.4.0] (2026-07-10)

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

## [0.3.0] (2026-05-14)

* Header-declared custom keys (`# CustomKeyDef=`) now drive typed parsing of entry values: registered keys are parsed via their `RegExp` (with pipe-split as fallback) and coerced per `FieldTypes` (`integer`, `decimal`, `boolean`, `date`, `time`, `string`, `enumeration`).
* Parsed custom-key values are exposed on `SequenceEntry.custom_values`, with the raw string preserved alongside for lossless round-trip writing.
* Undeclared custom keys continue to land in `SequenceEntry.extra` unchanged.
* `CustomKeyDef` now captures the `ConceptCURIE` field.
* Multiple `# CustomKeyDef=` lines per database are preserved instead of silently overwriting each other.
* Field tokenization is now quote-aware so a `RegExp` value containing escaped pipes survives splitting.

## [0.2.0] (2026-03-18)

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

## [0.1.0] (2026-03-18)

* First release on PyPI.
