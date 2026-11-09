"""Execute every docstring example in the package.

A ``>>>`` example is the only documentation the suite can prove is still true, so
every one of them runs here. The walk is deliberate: naming modules by hand is how
examples end up watched by nothing, and an import failure must break this test
rather than drop a module out of the count.
"""

from __future__ import annotations

import doctest
import importlib
import pkgutil
from types import ModuleType

import pefftacular

# Measured 2026-11-08 on Python 3.13.7: the walk finds 9 modules and they hold 0
# ``>>>`` examples between them. Raise these when you add modules or examples; never
# lower them without saying why in the commit message.
#
# MIN_EXAMPLES is 0 and that assert is inert on purpose - it arms itself the moment
# the first example lands. MIN_MODULES is the floor doing the work today: a walk that
# silently collects nothing also reports zero failures, and only this count notices.
MIN_MODULES = 9
MIN_EXAMPLES = 0


def _all_modules() -> list[ModuleType]:
    """Import every module in the package, the package itself included.

    No module is excluded by name: pefftacular has no runtime dependencies, so every
    module imports in any environment the tests run in. A module that needed an
    optional extra would be listed here explicitly, with a comment, so the exclusion
    stays visible.
    """
    mods = [pefftacular]
    for _finder, name, _ispkg in pkgutil.walk_packages(pefftacular.__path__, prefix=f"{pefftacular.__name__}."):
        mods.append(importlib.import_module(name))
    return mods


def test_docstring_examples() -> None:
    mods = _all_modules()
    assert len(mods) >= MIN_MODULES, f"module walk found only {len(mods)}"
    attempted = failed = 0
    for mod in mods:
        result = doctest.testmod(mod, verbose=False)
        attempted += result.attempted
        failed += result.failed
    assert failed == 0, f"{failed} of {attempted} docstring examples failed"
    assert attempted >= MIN_EXAMPLES, f"only {attempted} examples ran, expected >= {MIN_EXAMPLES}"
