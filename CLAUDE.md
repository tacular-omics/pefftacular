# CLAUDE.md

This project's guidance for AI coding agents lives in **[AGENTS.md](AGENTS.md)** —
read it first. It covers the module map, the `just check` commit gate, the
error/warning/logging model, and the format gotchas.

Quick reminders for Claude Code specifically:

- Run `just check` before committing; use `just fix` to auto-format.
- The PEFF spec is in `PEFF_SpecDoc_1.0_FINAL.pdf`; cite sections in comments.
- Parsing is permissive: warn via `PeffWarning`, never raise, on recoverable
  spec violations. Every new raise should carry a `hint=`.
