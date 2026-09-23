# Copilot instructions

The canonical guide for this repo is [`CLAUDE.md`](../CLAUDE.md) (commands, module map,
conventions, gotchas). Usage reference: [`llms-full.txt`](../llms-full.txt). Read those first.

Hard rules:

1. Use `just` recipes (`just --list`); fall back to `uv run ...`. Never install with pip,
   npm or another package manager. `just check` (format check, ruff, ty, pytest) must pass
   before any commit.
2. Keep the library dependency-free: `dependencies = []` in `pyproject.toml`.
3. Parsing is permissive: a recoverable spec violation is a
   `warnings.warn(msg, PeffWarning)`, never a raise. Every new raise is
   `PeffParseError`/`PeffWriteError` with a `hint=`.
4. Models are frozen, slotted dataclasses; do not mutate them. Import from the package
   root, not from `_`-prefixed modules.
5. Follow `PEFF_SpecDoc_1.0_FINAL.pdf` and cite the spec section in comments for new
   user-facing behaviour. Do not bump versions, tag or publish.
