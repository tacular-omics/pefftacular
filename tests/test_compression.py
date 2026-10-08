"""The ``compression`` argument of read_peff, PeffReader, to_records and write_peff."""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import typing
import warnings
from collections.abc import Callable
from pathlib import Path

import pytest

import pefftacular
from pefftacular import (
    Compression,
    FileHeader,
    PeffError,
    PeffParseError,
    PeffReader,
    SequenceEntry,
    read_peff,
    to_records,
    write_peff,
)

FIXTURES = Path(__file__).parent / "fixtures"
TEXT = (FIXTURES / "complex.peff").read_text(encoding="utf-8")
DATA = TEXT.encode()

COMPRESS: dict[str, Callable[[bytes], bytes]] = {"gzip": gzip.compress, "bz2": bz2.compress, "xz": lzma.compress}
DECOMPRESS: dict[str, Callable[[bytes], bytes]] = {
    "gzip": gzip.decompress,
    "bz2": bz2.decompress,
    "xz": lzma.decompress,
}
SUFFIX = {"gzip": ".gz", "bz2": ".bz2", "xz": ".xz"}
MAGIC = {"gzip": b"\x1f\x8b", "bz2": b"BZh", "xz": b"\xfd7zXZ\x00"}
KINDS = list(COMPRESS)

Parsed = tuple[FileHeader, list[SequenceEntry]]


@pytest.fixture(autouse=True)
def _quiet() -> typing.Iterator[None]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


@pytest.fixture(scope="module")
def plain() -> Parsed:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_peff(FIXTURES / "complex.peff")


def _kind_of(data: bytes) -> str | None:
    return next((k for k, m in MAGIC.items() if data.startswith(m)), None)


def test_compression_alias_is_exported() -> None:
    assert "Compression" in pefftacular.__all__
    assert set(typing.get_args(typing.get_args(Compression)[0])) == {"infer", "gzip", "bz2", "xz"}
    assert type(None) in typing.get_args(Compression)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("name", ["x.peff{}", "no_suffix{}.peff", "wrong.peff.txt"])
def test_read_infer_sniffs_magic_bytes(tmp_path: Path, plain: Parsed, kind: str, name: str) -> None:
    path = tmp_path / name.format(SUFFIX[kind])
    path.write_bytes(COMPRESS[kind](DATA))
    assert read_peff(path) == plain
    assert read_peff(path, compression="infer") == plain


@pytest.mark.parametrize("kind", KINDS)
def test_read_infer_plain_with_compressed_suffix(tmp_path: Path, plain: Parsed, kind: str) -> None:
    path = tmp_path / f"plain.peff{SUFFIX[kind]}"
    path.write_bytes(DATA)
    assert read_peff(path, compression="infer") == plain


@pytest.mark.parametrize("kind", KINDS)
def test_read_explicit_ignores_suffix(tmp_path: Path, plain: Parsed, kind: str) -> None:
    path = tmp_path / "x.peff"
    path.write_bytes(COMPRESS[kind](DATA))
    assert read_peff(path, compression=kind) == plain  # ty: ignore[invalid-argument-type]
    with PeffReader(path, compression=kind) as reader:  # ty: ignore[invalid-argument-type]
        assert (reader.header, list(reader)) == plain


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("actual", [None, *KINDS])
def test_read_explicit_rejects_other_format(tmp_path: Path, kind: str, actual: str | None) -> None:
    if actual == kind:
        return
    path = tmp_path / f"x.peff{SUFFIX[kind]}"
    path.write_bytes(DATA if actual is None else COMPRESS[actual](DATA))
    with pytest.raises(PeffParseError, match=f"not {kind}-compressed"):
        read_peff(path, compression=kind)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize("kind", KINDS)
def test_read_none_reads_compressed_bytes_as_plain(tmp_path: Path, kind: str) -> None:
    # None never decompresses, so compressed bytes fail as undecodable or as non-PEFF text.
    path = tmp_path / f"x.peff{SUFFIX[kind]}"
    path.write_bytes(COMPRESS[kind](DATA))
    with pytest.raises(PeffParseError):
        read_peff(path, compression=None)


def test_read_none_plain_file_with_gz_suffix(tmp_path: Path, plain: Parsed) -> None:
    path = tmp_path / "x.peff.gz"
    path.write_bytes(DATA)
    assert read_peff(path, compression=None) == plain


def test_read_text_handle_infer_and_none(plain: Parsed) -> None:
    assert read_peff(io.StringIO(TEXT)) == plain
    assert read_peff(io.StringIO(TEXT), compression=None) == plain


@pytest.mark.parametrize("kind", KINDS)
def test_read_text_handle_explicit_raises(kind: str) -> None:
    with pytest.raises(PeffError, match="needs a binary handle") as info:
        read_peff(io.StringIO(TEXT), compression=kind)  # ty: ignore[invalid-argument-type]
    assert any('"rb"' in note for note in info.value.__notes__)


@pytest.mark.parametrize("compression", ["infer", None])
def test_read_binary_handle_plain(plain: Parsed, compression: Compression) -> None:
    handle = io.BytesIO(DATA)
    assert read_peff(handle, compression=compression) == plain
    assert not handle.closed


@pytest.mark.parametrize("kind", KINDS)
def test_read_binary_handle_infer_does_not_decompress(kind: str) -> None:
    # "infer" on a handle reads it as it is, like a text handle.
    with pytest.raises(PeffParseError):
        read_peff(io.BytesIO(COMPRESS[kind](DATA)))


@pytest.mark.parametrize("kind", KINDS)
def test_read_binary_handle_explicit(plain: Parsed, kind: str) -> None:
    handle = io.BytesIO(COMPRESS[kind](DATA))
    assert read_peff(handle, compression=kind) == plain  # ty: ignore[invalid-argument-type]
    assert not handle.closed


@pytest.mark.parametrize("kind", KINDS)
def test_read_binary_handle_explicit_wrong_format(kind: str) -> None:
    with pytest.raises(PeffParseError, match=f"not {kind}-compressed"):
        read_peff(io.BytesIO(DATA), compression=kind)  # ty: ignore[invalid-argument-type]


@pytest.mark.parametrize("kind", KINDS)
def test_to_records_compression(tmp_path: Path, plain: Parsed, kind: str) -> None:
    path = tmp_path / "x.peff"
    path.write_bytes(COMPRESS[kind](DATA))
    expected = [e.to_record() for e in plain[1]]
    assert [r["sequence"] for r in to_records(path, compression=kind)] == [r["sequence"] for r in expected]  # ty: ignore[invalid-argument-type]
    with pytest.raises(PeffParseError):
        to_records(path, compression=None)


@pytest.mark.parametrize("bad", ["gz", "GZIP", "zip", "", 1])
def test_read_unknown_compression(tmp_path: Path, bad: object) -> None:
    path = tmp_path / "x.peff"
    path.write_bytes(DATA)
    for call in (
        lambda: read_peff(path, compression=bad),  # ty: ignore[invalid-argument-type]
        lambda: PeffReader(path, compression=bad),  # ty: ignore[invalid-argument-type]
        lambda: to_records(path, compression=bad),  # ty: ignore[invalid-argument-type]
    ):
        with pytest.raises(PeffError, match="valid values are 'infer', 'gzip', 'bz2', 'xz', None"):
            call()


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("upper", [False, True])
def test_write_infer_by_suffix(tmp_path: Path, plain: Parsed, kind: str, upper: bool) -> None:
    suffix = SUFFIX[kind].upper() if upper else SUFFIX[kind]
    path = tmp_path / f"x.peff{suffix}"
    write_peff(*plain, path)
    data = path.read_bytes()
    assert _kind_of(data) == kind
    assert read_peff(path) == plain


def test_write_infer_plain_suffix(tmp_path: Path, plain: Parsed) -> None:
    path = tmp_path / "x.peff"
    write_peff(*plain, path, compression="infer")
    assert path.read_bytes().startswith(b"# PEFF")


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("name", ["x.peff", "x.peff.gz", "x.peff.bz2", "x.peff.xz"])
def test_write_explicit_overrides_suffix(tmp_path: Path, plain: Parsed, kind: str, name: str) -> None:
    path = tmp_path / name
    write_peff(*plain, path, compression=kind)  # ty: ignore[invalid-argument-type]
    assert _kind_of(path.read_bytes()) == kind
    assert read_peff(path) == plain


@pytest.mark.parametrize("name", ["x.peff.gz", "x.peff.bz2", "x.peff.xz", "x.peff"])
def test_write_none_is_plain(tmp_path: Path, plain: Parsed, name: str) -> None:
    path = tmp_path / name
    write_peff(*plain, path, compression=None)
    assert path.read_bytes().startswith(b"# PEFF")
    assert read_peff(path) == plain


@pytest.mark.parametrize("compression", ["gzip", "infer"])
def test_write_gzip_is_reproducible(tmp_path: Path, plain: Parsed, compression: Compression) -> None:
    a, b = tmp_path / "a.peff.gz", tmp_path / "a2.peff.gz"
    write_peff(*plain, a, compression=compression)
    write_peff(*plain, b, compression=compression)
    # Same bytes apart from the file name stored in the gzip header.
    assert gzip.decompress(a.read_bytes()) == gzip.decompress(b.read_bytes())
    assert a.read_bytes()[4:8] == b"\x00\x00\x00\x00"  # mtime=0
    buf1, buf2 = io.BytesIO(), io.BytesIO()
    write_peff(*plain, buf1, compression="gzip")
    write_peff(*plain, buf2, compression="gzip")
    assert buf1.getvalue() == buf2.getvalue()


@pytest.mark.parametrize("compression", ["infer", None])
def test_write_text_handle_plain(plain: Parsed, compression: Compression) -> None:
    buf = io.StringIO()
    write_peff(*plain, buf, compression=compression)
    assert read_peff(io.StringIO(buf.getvalue())) == plain


@pytest.mark.parametrize("kind", KINDS)
def test_write_text_handle_explicit_raises(plain: Parsed, kind: str) -> None:
    buf = io.StringIO()
    with pytest.raises(PeffError, match="needs a binary handle") as info:
        write_peff(*plain, buf, compression=kind)  # ty: ignore[invalid-argument-type]
    assert type(info.value) is PeffError
    assert any("open in binary mode" in note and '"wb"' in note for note in info.value.__notes__)
    assert buf.getvalue() == ""


@pytest.mark.parametrize("compression", ["infer", None])
def test_write_binary_handle_plain(plain: Parsed, compression: Compression) -> None:
    buf = io.BytesIO()
    write_peff(*plain, buf, compression=compression)
    assert not buf.closed
    assert buf.getvalue().startswith(b"# PEFF")
    assert read_peff(io.BytesIO(buf.getvalue())) == plain


@pytest.mark.parametrize("kind", KINDS)
def test_write_binary_handle_explicit(tmp_path: Path, plain: Parsed, kind: str) -> None:
    buf = io.BytesIO()
    write_peff(*plain, buf, compression=kind)  # ty: ignore[invalid-argument-type]
    assert not buf.closed
    data = buf.getvalue()
    assert _kind_of(data) == kind
    assert DECOMPRESS[kind](data).decode().startswith("# PEFF")
    assert read_peff(io.BytesIO(data), compression=kind) == plain  # ty: ignore[invalid-argument-type]
    # A real file opened "wb" works the same and stays open.
    path = tmp_path / "handle.peff"
    with path.open("wb") as fh:
        write_peff(*plain, fh, compression=kind)  # ty: ignore[invalid-argument-type]
        assert not fh.closed
    assert read_peff(path) == plain


@pytest.mark.parametrize("bad", ["gz", "infer ", "zip", 0])
def test_write_unknown_compression(tmp_path: Path, plain: Parsed, bad: object) -> None:
    path = tmp_path / "x.peff"
    with pytest.raises(PeffError, match="valid values are 'infer', 'gzip', 'bz2', 'xz', None"):
        write_peff(*plain, path, compression=bad)  # ty: ignore[invalid-argument-type]
    assert not path.exists()


@pytest.mark.parametrize("compression", [None, "gzip", "bz2", "xz"])
def test_failed_write_leaves_binary_handle_open_and_untouched(plain: Parsed, compression: Compression) -> None:
    import dataclasses
    import gc

    bad = dataclasses.replace(plain[1][0], pname="foo\ud800")
    buf = io.BytesIO()
    with pytest.raises(UnicodeEncodeError):
        write_peff(plain[0], [bad], buf, compression=compression, verify=False)
    gc.collect()
    assert not buf.closed
    assert buf.getvalue() == b""
    write_peff(*plain, buf, compression=compression)  # still usable
    assert read_peff(io.BytesIO(buf.getvalue()), compression=compression) == plain


@pytest.mark.parametrize("kind", KINDS)
def test_infer_on_compressed_binary_handle_hints_explicit(kind: str) -> None:
    with pytest.raises(PeffParseError) as info:
        read_peff(io.BytesIO(COMPRESS[kind](DATA)))
    assert any("compression='gzip' (or bz2/xz)" in note for note in info.value.__notes__)
