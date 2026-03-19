"""Smoke test: verify public API is importable."""

from pefftacular import (
    FileHeader,
    PeffReader,
    SequenceEntry,
    read_peff,
    write_peff,
)


def test_public_api_importable():
    assert FileHeader is not None
    assert PeffReader is not None
    assert SequenceEntry is not None
    assert read_peff is not None
    assert write_peff is not None
