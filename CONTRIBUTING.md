# Contributing to pefftacular

Bug reports, fixes and new features are welcome. Use the
[issue tracker](https://github.com/tacular-omics/pefftacular/issues) for bugs and
questions, with a minimal example and the pefftacular and Python versions. Report
security problems privately as described in [SECURITY.md](SECURITY.md).

## Claiming an issue

Comment on the issue to say you are taking it before you start, and check that nobody has
already done so. One issue, one contributor: when two pull requests arrive for the same issue,
the first one that was claimed is reviewed and the other is closed, however good it is. If a
claim has had no pull request or update for 14 days, comment to take it over.

Issues labelled `good first issue` need no knowledge of the codebase. Issues labelled
`help wanted` do: each one names the files to read first and the test that shows it is done.

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

`just install` installs the whole test suite. pefftacular declares **no optional extras**
(its runtime dependencies are empty, so there is nothing to add to `uv sync`), and the
optional parts of the suite are test-only cross-checks against `fastatacular`,
`peptacular`, `pandas` and `polars` that live in the `dev` dependency group, which a bare
`uv sync` installs. The `docs` group is separate and no test needs it.

If you run the tests in some other environment, `uv run pytest tests` prints the answer
in its own header: `collected N items / M skipped` means M files were skipped at
collection, and a header with no `skipped` clause means none were. Under `-q` that clause
is not printed, so fall back to a subtraction: sum `passed`, `failed`, `skipped`,
`xfailed`, `xpassed` and `errors` from the summary line, leave out any subtest count, and
subtract the collected total. Do not read the skip count on its own, because both kinds of
guard move it: a `pytest.importorskip` inside a test function turns that one test into a
skip and leaves the collected total alone, while one at module import reports a single
skip standing for the whole file, whose tests are never collected at all.

Neither the header nor the subtraction sees a file dropped by a collect-time filter such
as `collect_ignore` or `--ignore`. For that, compare the test files on disk against the
files that produced test ids:

```bash
comm -23 <(ls tests/test_*.py | sort) \
  <(uv run pytest tests --co -q | grep '::' | sed 's/::.*//' | sort -u)
```

Anything it prints collected nothing, which is either a module skipped at import or a file
with no test functions in it. Every guard in this suite is inside a test function today, so
all three checks should come back clean. `uv pip list` tells you which of the four packages
you actually have.

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
