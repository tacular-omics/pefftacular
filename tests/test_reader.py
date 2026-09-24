"""Tests for the PeffReader resource contract (same shape as fastatacular.FastaReader)."""

from io import StringIO
from pathlib import Path

import pytest

from pefftacular import (
    CustomKeyValue,
    DatabaseHeader,
    FileHeader,
    PeffReader,
    PeffWarning,
    SequenceEntry,
    read_peff,
)

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


class TestReenter:
    def test_path_reader_restarts_on_second_with(self) -> None:
        # The header parsed in the first block must not pin the closed file handle.
        reader = PeffReader(FIXTURES / "minimal.peff")
        with reader:
            first_header = reader.header
        with reader:
            assert reader.header == first_header
            assert len(list(reader)) == 2

    def test_stream_reader_continues_on_second_with(self) -> None:
        # A stream is not reopened: the second block continues where the first stopped.
        reader = PeffReader(StringIO((FIXTURES / "minimal.peff").read_text()))
        with reader:
            header = reader.header
        with reader:
            assert reader.header == header
            assert len(list(reader)) == 2


_COMMENT_BODY = """\
# PEFF 1.0
# //
# DbName=t
# Prefix=sp
# DbVersion=1
# DbSource=x
# NumberOfEntries=2
# SequenceType=AA
# //
>sp:P1 \\PName=one
; a FASTA comment line
ACDE
;another
FGHI
>sp:P2
KLMN
"""


class TestCommentLines:
    def test_semicolon_lines_skipped(self) -> None:
        _, entries = read_peff(StringIO(_COMMENT_BODY))
        assert [e.sequence for e in entries] == ["ACDEFGHI", "KLMN"]

    def test_semicolon_line_before_first_entry_skipped(self) -> None:
        text = _COMMENT_BODY.replace(">sp:P1", "; leading comment\n>sp:P1", 1)
        _, entries = read_peff(StringIO(text))
        assert entries[0].sequence == "ACDEFGHI"

    def test_hash_line_in_entry_body_skipped_with_warning(self) -> None:
        text = _COMMENT_BODY.replace(";another\n", "# DbName=late\n")
        with pytest.warns(PeffWarning, match="line 13"):
            _, entries = read_peff(StringIO(text))
        assert entries[0].sequence == "ACDEFGHI"


class TestHeaderModelsHash:
    @pytest.mark.parametrize(
        "obj",
        [
            DatabaseHeader(prefix="sp"),
            FileHeader(peff_version="1.0"),
            CustomKeyValue(key_name="K", fields={"a": 1}),
        ],
    )
    def test_unhashable_but_comparable(self, obj: object) -> None:
        assert type(obj).__hash__ is None
        assert obj == obj
        with pytest.raises(TypeError):
            hash(obj)
