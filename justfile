# pefftacular task runner.
#
# For humans and coding agents alike. The one command to run before committing
# is `just check` — it runs the exact gate CI enforces (format, lint, types,
# tests) read-only. Use `just fix` to auto-apply formatting/lint fixes.
#
# Run `just` (or `just --list`) to see every recipe.

# List all recipes (default)
default:
    @just --list

# Install project + dev dependencies into the uv-managed venv
install:
    uv sync

# Alias for install
dev:
    uv sync

# --- The commit gate ------------------------------------------------------

# Full read-only gate: formatting, lint, types, and tests (run this before committing)
check:
    uv run ruff format --check src tests
    uv run ruff check src tests
    uv run ty check src
    uv run pytest tests

# Auto-fix imports, lint, and formatting across src and tests
fix:
    uv run ruff check --fix src tests
    uv run ruff format src tests

# --- Individual steps -----------------------------------------------------

# Lint (no changes) across src and tests
lint:
    uv run ruff check src tests

# Format code across src and tests
format:
    uv run ruff check --select I --fix src tests
    uv run ruff format src tests

# Type-check the package with ty
ty:
    uv run ty check src

# Run the test suite
test:
    uv run pytest tests

# Run the test suite with verbose output
test-v:
    uv run pytest tests -v

# Run tests for a specific file, e.g. `just test-file tests/test_errors.py`
test-file FILE:
    uv run pytest {{FILE}} -v

# Run tests with coverage (terminal report)
cov:
    uv run pytest tests --cov=src/pefftacular --cov-report=term-missing

# Run tests with coverage (XML for Codecov)
test-cov:
    uv run pytest tests --cov=src/pefftacular --cov-report=xml

# Run tests with JUnit XML report for Codecov test-results upload
codecov-tests:
    uv run pytest tests --junit-xml=junit.xml

# --- Packaging & housekeeping --------------------------------------------

# Build the sdist + wheel
build:
    uv build

# Remove caches and compiled files
clean:
    find . -type d -name __pycache__ -exec rm -rf {} +
    find . -type d -name .pytest_cache -exec rm -rf {} +
    find . -type d -name .ruff_cache -exec rm -rf {} +
    find . -name "*.pyc" -delete

# Serve docs locally (requires mkdocs)
docs:
    uv run mkdocs serve

# Deploy docs to GitHub Pages (requires mkdocs)
docs-deploy:
    uv run mkdocs gh-deploy --force
