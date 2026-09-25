"""Shared pytest configuration: Hypothesis profiles.

Pick one with ``HYPOTHESIS_PROFILE``:

- ``default`` (30 examples): the local run, a few seconds for every property test.
- ``ci`` (100 examples, 3000 for the lexer fast-path checks, 500 for the escaped
  round trip): what CI runs, the counts the suite used before the default was cut.
- ``thorough`` (5000 examples): run it locally before touching the lexer, parser
  or writer:

    HYPOTHESIS_PROFILE=thorough uv run pytest tests/test_properties.py

Tests that need more examples than the profile keep their higher count through
``settings(max_examples=...)`` only when the profile is not ``default``.
"""

import os

from hypothesis import HealthCheck, Phase, settings

# The ``explain`` phase traces every line of the parser and takes minutes per failure.
_PHASES = tuple(p for p in Phase if p is not Phase.explain)

settings.register_profile("default", max_examples=30, deadline=None, phases=_PHASES)
settings.register_profile("ci", max_examples=100, deadline=None, phases=_PHASES)
settings.register_profile(
    "thorough",
    max_examples=5_000,
    deadline=None,
    phases=_PHASES,
    suppress_health_check=[HealthCheck.too_slow],
)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
