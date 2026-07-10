# AGENTS.md

Guidance for AI coding agents (and humans) working in this repository. If you
change how the project is built, tested, or structured, update this file.

## What this is

`pefftacular` is a **pure-Python, zero-runtime-dependency** library for reading
and writing [PEFF](http://www.psidev.info/peff) (PSI Extended FASTA Format)
files — FASTA plus rich per-entry proteomics annotations (PTMs, variants,
processed forms, proteoforms, disulfide bonds) encoded in the description line.

The authoritative format spec is in the repo: **`PEFF_SpecDoc_1.0_FINAL.pdf`**.
When in doubt about parsing/writing behavior, read the relevant section (the
code comments cite section numbers, e.g. "spec §3.3.3").

- Python **3.12+** only. Target `py312`.
- Lint/format: **ruff** (line length **120**). Type-check: **ty**.
- Package layout: `src/` (src layout), managed with **uv**.

## Commands

Run these through `just` (see `justfile` for the full list):

| Task | Command |
| --- | --- |
| Install deps | `just install` |
| **Pre-commit gate (do this before committing)** | `just check` |
| Auto-fix lint + formatting | `just fix` |
| Tests only | `just test` |
| One test file | `just test-file tests/test_errors.py` |
| Coverage | `just cov` |

`just check` runs, read-only and in this order: `ruff format --check` →
`ruff check` → `ty check src` → `pytest`. It must pass before any commit.

## Module map (`src/pefftacular/`)

| Module | Responsibility |
| --- | --- |
| `_models.py` | Frozen `@dataclass(slots=True)` models: `FileHeader`, `DatabaseHeader`, `SequenceEntry`, and the annotation types (`VariantSimple`, `ModResUnimod`, `DisulfideBond`, `Proteoform`, …). Pure data, no logic. |
| `_lexer.py` | Depth-/escape-aware tokenizers for the description line: `split_items` (parenthesized items), `split_fields` (pipe-separated components), `split_description_keys` (the `\Key=value` scan). |
| `_parser.py` | `PeffReader` (lazy) + `read_peff` (eager). Header parsing, per-key annotation parsing, and spec validation. |
| `_writer.py` | `write_peff` — serializes models back to canonical PEFF text. |
| `errors.py` | `PeffError` base, `PeffParseError`, `PeffWriteError`, and the `PeffWarning` category. |
| `__init__.py` | Public API surface (`__all__`) + library `NullHandler` for logging. |

The public API is exactly what `__init__.py` exports. Underscore-prefixed
modules are internal; import from the package root in tests and examples.

## How the library talks to you (errors, warnings, logging)

This is the machinery to lean on when debugging — it is designed to make the
library's behavior legible.

**Exceptions.** Everything derives from `PeffError` (itself a `ValueError`).
`PeffParseError` carries structured, machine-readable attributes:
- `.line` — 1-based absolute line number (when known)
- `.context` — the exact offending text
- `.hint` — a short, actionable fix

`.context` and `.hint` are also attached as PEP-678 exception **notes**, so they
appear in tracebacks without polluting `str(err)`. `PeffWriteError` carries
`.hint` too. Prefer catching `PeffError` for "any pefftacular failure".

**Warnings.** Parsing is **deliberately permissive** — it returns the data even
when the file violates the spec — but every spec `MUST`-violation is emitted
through the `PeffWarning` category (a `UserWarning` subclass): out-of-range
positions, missing required fields, header/entry-count mismatches, un-coercible
custom values, etc. To make these fatal:

```python
import warnings
from pefftacular import PeffWarning
warnings.simplefilter("error", PeffWarning)
```

**Logging.** Standard library-logging pattern — the package attaches a
`NullHandler` and never configures logging itself. Turn on a behavioral trace:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("pefftacular").setLevel(logging.DEBUG)
```

Loggers: `pefftacular.parser` and `pefftacular.writer`. `INFO` marks milestones
(entries read/written); `DEBUG` traces file open, header parse, and entry counts.

## Conventions

- **New raises:** always pass a `hint=` (and `context=`/`line=` where relevant).
- **New warnings:** use `warnings.warn(msg, PeffWarning, stacklevel=...)`. Never
  raise for a recoverable spec violation — warn and keep the data.
- **New user-facing behavior:** cite the spec section in a comment.
- Models are frozen and slotted — don't mutate; construct new instances.
- Keep the library dependency-free (`dependencies = []` in `pyproject.toml`).

## Gotchas (learned the hard way — don't regress these)

- **Escaping (spec §3.3.3).** In entry items, `\|`, `\(`, `\)`, `\\` are literal;
  *balanced* parens (e.g. `N-linked (GlcNAc...)`) are **not** escaped. The
  writer re-escapes on output. Entry-item tokenizing uses
  `split_fields(item, unescape=True)`; `CustomKeyDef` header values use the
  default quoted mode (`unescape=False`) so regex backslashes survive.
- **`DisulfideBond.annot_id_refs`** holds annotation-ID references to prior
  `ModResPsi` entries — **not** residue positions (spec §3.4.2). That's why it's
  excluded from position-range validation.
- **Header flag casing.** `ProteoformDb`/`HasAnnotationIdentifiers` are parsed
  case-insensitively because the spec's own examples disagree with its text.
- **Round-trip tests** compare parsed *models*, not raw text — escape/unescape
  must be exact inverses.

## Testing

- Tests live in `tests/`, one file per area; fixtures in `tests/fixtures/`.
- Assert warnings with `pytest.warns(PeffWarning, match=...)`; assert *absence*
  with `warnings.simplefilter("error", PeffWarning)` inside `catch_warnings`.
- A "valid" fixture must parse **without** emitting a `PeffWarning`. If you add
  annotations to a fixture, keep positions within `1..len(sequence)`.
- Add or update tests for every behavior change; keep `just check` green.

## Releasing

Version lives in `pyproject.toml`; changelog in `HISTORY.md`. Pre-1.0, breaking
changes go in a minor bump (e.g. 0.3 → 0.4). Work on a `release/x.y.z` branch,
finalize the `HISTORY.md` heading with the date, and keep `just check` green.
