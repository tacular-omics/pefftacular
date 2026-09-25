# Contributing to pefftacular

Bug reports, fixes and new features are welcome. Use the
[issue tracker](https://github.com/tacular-omics/pefftacular/issues) for bugs and
questions, with a minimal example and the pefftacular and Python versions. Report
security problems privately as described in [SECURITY.md](SECURITY.md).

## Development setup

You need Python 3.12 or newer, [uv](https://docs.astral.sh/uv/) and
[just](https://just.systems/).

```bash
git clone https://github.com/tacular-omics/pefftacular.git
cd pefftacular
just install
```

`just --list` shows the other recipes.

Install the commit hooks once per clone with `uvx pre-commit install`. On every
commit they run `ruff check`, `ruff format --check` and a few file checks
(`uvx pre-commit run --all-files` runs them on everything). Type checks and
tests are not in the hooks; run `just check` for those.

## Running checks

CI runs these on every pull request. Run them before you push:

```bash
uv run ruff check src tests
uv run ruff format --check src tests
uv run ty check src
uv run pytest tests
python scripts/release_version.py check
```

`uv run pytest tests` uses a small Hypothesis profile (30 examples) so it finishes in
seconds. CI runs `HYPOTHESIS_PROFILE=ci` (100 examples, more for the lexer checks);
`HYPOTHESIS_PROFILE=thorough` (5000) is for changes to the lexer, parser or writer.

If you change the documentation, build it the way CI does:

```bash
uv run mkdocs build --strict
```

## Pull requests

- Branch from `main` and keep each pull request to one change.
- Add tests for new behaviour and fixed bugs.
- Add a line under `## [Unreleased]` in `CHANGELOG.md` for user-visible changes.
- Do not change the version number; the maintainers set it when releasing.
- CI must pass: lint, type check, tests on Python 3.12 to 3.14, and the wheel build.

## CLAUDE.md and AI agents

`CLAUDE.md` is the maintained developer guide: commands, architecture,
conventions and known gotchas. `AGENTS.md` points other coding agents to it.
Read it before a larger change, and update it in the same pull request when your
change makes part of it wrong. AI-assisted contributions are welcome; you are
responsible for reviewing and testing what you submit.
