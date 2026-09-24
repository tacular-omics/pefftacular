"""Tests for the PeffReader resource contract (same shape as fastatacular.FastaReader)."""

from io import StringIO
from pathlib import Path

import pytest

from pefftacular import PeffReader, SequenceEntry, read_peff

FIXTURES = Path(__file__).parent / "fixtures"


class TestContextManager:
    def test_init_does_not_open_file(self, tmp_path: Path) -> None:
        # Constructing a reader for a missing path must not touch the filesystem.
        reader = PeffReader(tmp_path / "missing.peff")
        with pytest.raises(FileNotFoundError):
            reader.__enter__()

    def test_no_file_handle_before_enter(self) -> None:
        reader = PeffReader(FIXTURES / "minimal.peff")
        assert reader._fh is None

    def test_iter_outside_with_raises(self) -> None:
        with pytest.raises(RuntimeError, match="context manager"):
            iter(PeffReader(FIXTURES / "minimal.peff"))

    def test_header_outside_with_raises(self) -> None:
        with pytest.raises(RuntimeError, match="context manager"):
            PeffReader(StringIO("# PEFF 1.0\n")).header  # noqa: B018

    def test_with_path(self) -> None:
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            assert reader.header.peff_version == "1.0"
            assert len(list(reader)) == 2

    def test_stream_not_closed(self) -> None:
        stream = StringIO((FIXTURES / "minimal.peff").read_text())
        with PeffReader(stream) as reader:
            list(reader)
        assert not stream.closed

    def test_owned_file_closed_on_exit(self) -> None:
        with PeffReader(FIXTURES / "minimal.peff") as reader:
            fh = reader._fh
            assert fh is not None
        assert fh.closed


class TestBom:
    def test_utf8_bom_path(self, tmp_path: Path) -> None:
        path = tmp_path / "bom.peff"
        path.write_bytes(b"\xef\xbb\xbf" + (FIXTURES / "minimal.peff").read_bytes())
        header, entries = read_peff(path)
        assert header.peff_version == "1.0"
        assert len(entries) == 2

    def test_bom_in_stream(self) -> None:
        text = "\ufeff" + (FIXTURES / "minimal.peff").read_text()
        header, _ = read_peff(StringIO(text))
        assert header.peff_version == "1.0"


class TestSequenceEntryHash:
    def test_unhashable_but_comparable(self) -> None:
        a = SequenceEntry(prefix="sp", db_unique_id="X", sequence="AC")
        b = SequenceEntry(prefix="sp", db_unique_id="X", sequence="AC")
        assert a == b
        assert SequenceEntry.__hash__ is None
        with pytest.raises(TypeError):
            hash(a)
