# pefftacular — Claude Code Guide

## Project overview

`pefftacular` is a **pure-Python, zero-runtime-dependency** library for reading and
writing [PEFF](https://www.psidev.info/peff) (PSI Extended FASTA Format) files: FASTA
plus rich per-entry proteomics annotations (PTMs, variants, processed forms,
proteoforms, disulfide bonds, header-declared custom keys) encoded in the description
line. Users are people building or consuming annotated protein sequence databases.

The authoritative format spec is in the repo: **`PEFF_SpecDoc_1.0_FINAL.pdf`**. When in
doubt about parsing/writing behaviour, read the relevant section; code comments cite
section numbers (e.g. "spec §3.3.3").

Place in the tacular-omics graph:

- **Upstream:** none. Tier 0, no sibling dependencies (`dependencies = []`).
- **Sibling:** `fastatacular` (plain FASTA) shares the same API shape
  (`read_*` eager, `*Reader` lazy, `write_*`).
- **Downstream (outside the core workspace):** `peff_digest`, `peff_uniprot_fetcher`.
  A breaking change to the models or reader must be noted for them.

Python **3.12+** (CI tests 3.12-3.14). Ruff line length **120**, target `py312`. Type
checker: **ty**. src layout, managed with **uv**, built with hatchling.

## Commands

```bash
just                 # list every recipe
just install         # uv sync (alias: just dev)
just check           # COMMIT GATE, read-only: ruff format --check, ruff check, ty check src, pytest
just fix             # ruff check --fix + ruff format (src, tests)
just lint            # ruff check src tests
just format          # ruff isort fix + ruff format
just ty              # ty check src
just test            # pytest tests
just test-v          # pytest tests -v
just test-file tests/test_errors.py   # one file, verbose
just cov             # pytest with --cov=src/pefftacular, term-missing report
just test-cov        # coverage XML (Codecov)
just build           # uv build (sdist + wheel)
just clean           # remove caches
just check-version   # fail if __version__ / CITATION.cff disagree
```

`just check` must pass before any commit (187 tests, well under a second).
`just docs` / `just docs-deploy` exist but are dead: there is no `mkdocs.yml` and mkdocs
is not a dependency. Do not rely on them.

## Architecture

```
src/pefftacular/
  __init__.py   # public API (__all__), __version__, library NullHandler for logging
  _models.py    # frozen @dataclass(slots=True) models; pure data, no logic
  _lexer.py     # depth-/escape-aware tokenizers for the description line:
                #   split_items (parenthesized items), split_fields (| components),
                #   split_description_keys (the "\Key=value" scan)
  _parser.py    # PeffReader (lazy) + read_peff (eager): header parsing, per-key
                #   annotation parsing, custom-key coercion, spec validation warnings
  _writer.py    # write_peff: serializes models back to canonical PEFF text
  errors.py     # PeffError base, PeffParseError, PeffWriteError, PeffWarning
scripts/release_version.py   # version sync/check used by the release recipes
tests/                       # one file per area; fixtures in tests/fixtures/*.peff
PEFF_SpecDoc_1.0_FINAL.pdf   # the spec
```

Data flow (read): `PeffReader(source)` → first access to `.header` or iteration runs
`_parse_file_header` (version line, `GeneralComment`s, one `DatabaseHeader` per `# //`
block) → each `>` line plus its sequence lines goes to `_parse_entry`, which calls
`split_description_keys`, dispatches each known key to its `_parse_*` helper (which uses
`split_items` + `split_fields(..., unescape=True)`), routes header-declared custom keys
to `_parse_custom_value`, and puts everything else in `extra` → `_warn_on_invalid_annotations`
→ `SequenceEntry`. After the last entry, `NumberOfEntries` is checked per prefix.

Data flow (write): `write_peff` validates entries (non-empty prefix / id / sequence),
writes the header, then each entry's `\Key=value` pairs in a fixed canonical order,
sequence wrapped at 60 characters.

## Public API

Exactly what `pefftacular.__all__` exports. Underscore modules are internal; import from
the package root in tests and examples.

- **I/O:** `read_peff(source)` → `(FileHeader, list[SequenceEntry])`;
  `PeffReader(source)` lazy reader (`.header`, iterate for entries, context manager);
  `write_peff(header, entries, dest)`.
- **Header models:** `FileHeader`, `DatabaseHeader`, `CustomKeyDef`, `OptionalTagDef`.
- **Entry model:** `SequenceEntry`.
- **Annotation models:** `VariantSimple`, `VariantComplex`, `ModResUnimod`, `ModResPsi`,
  `ModRes`, `Processed`, `DisulfideBond`, `Proteoform`, `SequenceRange`, `CustomKeyValue`.
- **Errors/warnings:** `PeffError` (base, subclasses `ValueError`), `PeffParseError`
  (`.line`, `.context`, `.hint`), `PeffWriteError` (`.hint`), `PeffWarning`
  (`UserWarning` subclass).
- `__version__`.

Full signatures and examples: `llms-full.txt`.

## Conventions

- **Docstrings:** Google style, short. Explain the why and cite the spec section for any
  user-facing behaviour.
- **Typing:** full hints, `X | None` not `Optional`, builtin generics. `ty check src`
  must be clean.
- **Models** are frozen and slotted: never mutate, construct new instances
  (`dataclasses.replace`).
- **Errors:** every new raise uses `PeffParseError` / `PeffWriteError` with a `hint=`
  (and `context=` / `line=` where known). `context` and `hint` are also attached as
  PEP 678 notes so they show in tracebacks without changing `str(err)`.
- **Warnings:** parsing is deliberately permissive. For a recoverable spec violation,
  `warnings.warn(msg, PeffWarning, stacklevel=...)` and keep the data; never raise.
- **Logging:** stdlib only. Loggers `pefftacular.parser` and `pefftacular.writer`;
  `INFO` for milestones (entries read/written), `DEBUG` for file open, header parse,
  counts. The package attaches a `NullHandler` and never configures logging itself.
- **Dependencies:** keep runtime dependencies empty.
- **Tests:** `tests/`, one file per area, fixtures in `tests/fixtures/`. Assert warnings
  with `pytest.warns(PeffWarning, match=...)`; assert their absence with
  `warnings.simplefilter("error", PeffWarning)` inside `warnings.catch_warnings()`.
  A "valid" fixture must parse without any `PeffWarning`; keep annotation positions
  within `1..len(sequence)`. Add or update tests for every behaviour change.

## Gotchas (learned the hard way; do not regress these)

- **Escaping (spec §3.3.3).** In entry items `\|`, `\(`, `\)`, `\\` are literal;
  *balanced* parens (e.g. `N-linked (GlcNAc...)`) are **not** escaped. The writer
  re-escapes on output (`_escape_component`, parens only when unbalanced). Entry items
  are tokenized with `split_fields(item, unescape=True)`; `CustomKeyDef` header values
  use the default quoted mode (`unescape=False`) so regex backslashes survive.
- **`split_description_keys` does not honour escapes.** It counts every `(` and `)`
  for depth, so an item containing an escaped unbalanced paren (`name \( foo`) leaves
  depth > 0 and swallows every following `\Key=` into that value. The writer produces
  exactly this for names with unbalanced parens. Known bug, not yet fixed.
- **`DisulfideBond.annot_id_refs`** holds annotation-ID references to prior
  `ModResPsi` entries, **not** residue positions (spec §3.4.2). That is why it is
  excluded from position-range validation.
- **Header flag casing.** `ProteoformDb` / `HasAnnotationIdentifiers` are parsed
  case-insensitively (also `IsProteoformDB`, `HasAnnotationIdentifier`) because the
  spec's own examples disagree with its text. The writer always emits
  `HasAnnotationIdentifiers=true` / `ProteoformDb=true`.
- **Round-trip tests** compare parsed *models*, not raw text: escape/unescape must be
  exact inverses, but key order, `\OX` → `\NcbiTaxId`, and header key order normalize.
- **`\Variant=` (deprecated)** emits `DeprecationWarning`, not `PeffWarning`, and the
  value lands in `extra["Variant"]`.
- **`CustomKeyValue.raw` wins on write.** If `raw` is non-empty the writer emits it
  verbatim and ignores `fields`; clear `raw` when building a value from edited fields.
- **Header keys `SpecificKey` / `SpecificValue`** are parsed and then dropped; they do
  not survive a round trip. Unknown single-valued header keys go to
  `DatabaseHeader.extra`.
- **Duplicate description keys:** the last `\Key=` on a line wins silently.
- **`PeffReader` iterates once.** It wraps a single line iterator; a second `for` over
  the same reader yields nothing. Use `with PeffReader(path)` so an owned file closes.
- **Free-text scalar values are not escaped by the writer** (`pname`, `gname`,
  `tax_name`, `comment`, `extra` values). A value containing ` \Key=` is read back as a
  new key.

## Releasing

Only the tacular-omics overseer bumps versions or publishes. See `just --list`
(`set-version`, `sync-version`, `check-version`) and the workspace release checklist.
Version source: `__version__` in `src/pefftacular/__init__.py` (hatch reads it), kept in
sync with `CITATION.cff` by `scripts/release_version.py`. Changelog: `CHANGELOG.md`
(`## [X.Y.Z] (YYYY-MM-DD)`, newest first, `[Unreleased]` on top). Pre-1.0, breaking
changes go in a minor bump. Publishing is the `release: published` GitHub workflow
(PyPI trusted publishing); Zenodo archives each GitHub release.

## Workspace note

This repo is also developed inside the tacular-omics uv workspace; there `uv run` uses
the shared `.venv`. See the workspace CLAUDE.md.
