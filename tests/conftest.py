"""Shared pytest configuration: Hypothesis profiles.

``default`` keeps the property tests to a few seconds. Run the ``thorough`` profile
locally before touching the lexer, parser or writer:

    HYPOTHESIS_PROFILE=thorough uv run pytest tests/test_properties.py
"""

import os

from hypothesis import HealthCheck, Phase, settings

# The ``explain`` phase traces every line of the parser and takes minutes per failure.
_PHASES = tuple(p for p in Phase if p is not Phase.explain)

settings.register_profile("default", max_examples=100, deadline=None, phases=_PHASES)
settings.register_profile(
    "thorough",
    max_examples=5_000,
    deadline=None,
    phases=_PHASES,
    suppress_health_check=[HealthCheck.too_slow],
)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
