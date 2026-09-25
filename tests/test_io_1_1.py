"""1.1 I/O: escape round trips, ``write_peff(verify=...)``, streaming writes, compressed input."""

from __future__ import annotations

import bz2
import dataclasses
import gzip
import io
import lzma
import os
import subprocess
import sys
import threading
import time
import warnings
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pefftacular import (
    DatabaseHeader,
    FileHeader,
    ModResPsi,
    PeffParseError,
    PeffReader,
    PeffWriteError,
    SequenceEntry,
    VariantSimple,
    read_peff,
    write_peff,
)
from pefftacular import _writer as writer

FIXTURES = Path(__file__).parent / "fixtures"
HEADER = FileHeader(peff_version="1.0", databases=(DatabaseHeader(prefix="sp"),))


def _write(entries, **kwargs) -> str:  # type: ignore[no-untyped-def]
    buf = io.StringIO()
    write_peff(HEADER, entries, buf, **kwargs)
    return buf.getvalue()


def _read(text: str) -> list[SequenceEntry]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_peff(io.StringIO(text))[1]


# ---------------------------------------------------------------------------
# Escape round trip: values full of PEFF syntax survive write -> read
# ---------------------------------------------------------------------------

_special = st.sampled_from(list("\\|()= ab") + ["\\\\", "\\(", "\\|", "(x|y)"])
_component = st.lists(_special, max_size=10).map("".join)
_nonempty = _component.filter(bool)


def _examples(n: int) -> settings:
    """``n`` examples under the ``ci``/``thorough`` profiles, the profile count locally."""
    if os.environ.get("HYPOTHESIS_PROFILE", "default") == "default":
        return settings()
    return settings(max_examples=max(n, settings().max_examples))


@given(tag=_nonempty, new_aa=_component, mod_name=_nonempty)
@_examples(500)
def test_escaped_components_round_trip(tag: str, new_aa: str, mod_name: str) -> None:
    entry = SequenceEntry(
        prefix="sp",
        db_unique_id="P1",
        sequence="MKTAYIAK",
        variant_simple=(VariantSimple(position=2, new_amino_acid=new_aa, tag=tag),),
        mod_res_psi=(ModResPsi(positions=(3,), accession="MOD:00046", name=mod_name),),
    )
    text = _write([entry])
    [back] = _read(text)
    assert back.variant_simple == entry.variant_simple
    assert back.mod_res_psi == entry.mod_res_psi
    # The unverified write of the same entry is byte-identical.
    assert _write([entry], verify=False) == text


# ---------------------------------------------------------------------------
# verify=False and streaming
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["complex.peff", "custom_keys.peff", "multidb.peff", "PEFF_AnnotID_Insulin_Valid.peff"]
)
def test_verify_false_writes_the_same_text(name: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        header, entries = read_peff(FIXTURES / name)
    a, b = io.StringIO(), io.StringIO()
    write_peff(header, entries, a)
    write_peff(header, entries, b, verify=False)
    assert a.getvalue() == b.getvalue()


def test_verify_true_still_catches_values_that_read_back_differently() -> None:
    bad = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK", extra={"Xk": "a \\GName=y"})
    with pytest.raises(PeffWriteError, match="gname does not read back"):
        _write([bad])
    # Skipping the check writes what the caller asked for; it reads back differently.
    assert _read(_write([bad], verify=False))[0].gname == "y"


def test_basic_checks_run_with_verify_false() -> None:
    with pytest.raises(PeffWriteError, match="empty prefix"):
        _write([SequenceEntry(prefix="", db_unique_id="P1", sequence="MK")], verify=False)


def test_entries_are_consumed_as_a_stream() -> None:
    seen: list[int] = []

    def gen() -> Iterator[SequenceEntry]:
        for i in range(5):
            seen.append(i)
            yield SequenceEntry(prefix="sp", db_unique_id=f"P{i}", sequence="MK")

    text = _write(gen())
    assert seen == [0, 1, 2, 3, 4]
    assert [e.db_unique_id for e in _read(text)] == [f"P{i}" for i in range(5)]


def test_error_mid_stream_writes_nothing(tmp_path: Path) -> None:
    def gen() -> Iterator[SequenceEntry]:
        yield SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK")
        yield SequenceEntry(prefix="", db_unique_id="P2", sequence="MK")

    path = tmp_path / "out.peff"
    with pytest.raises(PeffWriteError, match="Entry 1"):
        write_peff(HEADER, gen(), path)
    assert not path.exists()
    buf = io.StringIO()
    with pytest.raises(PeffWriteError):
        write_peff(HEADER, gen(), buf)
    assert buf.getvalue() == ""


def test_write_spills_to_disk_past_the_spool_limit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(writer, "_SPOOL_BYTES", 64)
    entries = [SequenceEntry(prefix="sp", db_unique_id=f"P{i}", sequence="MKTAYIAK" * 20) for i in range(50)]
    path = tmp_path / "big.peff"
    write_peff(HEADER, entries, path)
    monkeypatch.setattr(writer, "_SPOOL_BYTES", 32 * 1024 * 1024)
    assert path.read_text() == _write(entries)
    assert [dataclasses.astuple(e) for e in read_peff(path)[1]] == [
        dataclasses.astuple(e) for e in _read(_write(entries))
    ]


# ---------------------------------------------------------------------------
# Compressed and undecodable input
# ---------------------------------------------------------------------------

COMPRESSORS: dict[str, Callable[[bytes], bytes]] = {".gz": gzip.compress, ".bz2": bz2.compress, ".xz": lzma.compress}
TEXT = (FIXTURES / "complex.peff").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def plain() -> tuple[FileHeader, list[SequenceEntry]]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_peff(FIXTURES / "complex.peff")


@pytest.mark.parametrize("suffix", list(COMPRESSORS))
@pytest.mark.parametrize("name", ["x.peff{}", "no_suffix_{}.peff"])
def test_compressed_input_reads_like_plain(tmp_path: Path, plain, suffix: str, name: str) -> None:  # type: ignore[no-untyped-def]
    path = tmp_path / name.format(suffix)
    path.write_bytes(COMPRESSORS[suffix](TEXT.encode()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert read_peff(path) == plain
        assert read_peff(str(path)) == plain
        with PeffReader(path) as reader:
            assert (reader.header, list(reader)) == plain


def test_compressed_with_bom(tmp_path: Path, plain) -> None:  # type: ignore[no-untyped-def]
    path = tmp_path / "bom.peff.gz"
    path.write_bytes(gzip.compress(("\ufeff" + TEXT).encode()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert read_peff(path) == plain


@pytest.mark.parametrize("suffix", list(COMPRESSORS))
def test_truncated_compressed_input_raises_peff_parse_error(tmp_path: Path, suffix: str) -> None:
    data = COMPRESSORS[suffix]((TEXT * 20).encode())
    path = tmp_path / f"cut.peff{suffix}"
    path.write_bytes(data[: len(data) // 2])
    with pytest.raises(PeffParseError, match="Cannot read the input") as info, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        read_peff(path)
    assert info.value.__cause__ is not None


def test_gz_suffix_on_plain_text_reads_as_plain(tmp_path: Path, plain) -> None:  # type: ignore[no-untyped-def]
    # The magic bytes decide, not the file name.
    path = tmp_path / "plain.peff.gz"
    path.write_bytes(TEXT.encode())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert read_peff(path) == plain


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs os.mkfifo")
@pytest.mark.parametrize("suffix", ["", *COMPRESSORS])
def test_fifo_is_read_once(tmp_path: Path, plain, suffix: str) -> None:  # type: ignore[no-untyped-def]
    # A pipe can be read only once: the format sniff must not consume its first bytes.
    data = COMPRESSORS[suffix](TEXT.encode()) if suffix else TEXT.encode()
    fifo = tmp_path / "pipe.peff"
    os.mkfifo(fifo)

    def feed() -> None:
        with fifo.open("wb") as w:
            w.write(data)

    def read() -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result.append(read_peff(fifo))

    result: list = []  # type: ignore[type-arg]
    writer_thread = threading.Thread(target=feed, daemon=True)
    reader_thread = threading.Thread(target=read, daemon=True)
    writer_thread.start()
    reader_thread.start()
    reader_thread.join(timeout=10)
    if reader_thread.is_alive():  # a second open() of the FIFO blocks: give it EOF, then fail
        os.close(os.open(fifo, os.O_WRONLY | os.O_NONBLOCK))
        pytest.fail("reading the FIFO blocked (was it opened twice?)")
    writer_thread.join(timeout=10)
    assert result == [plain]


@pytest.mark.skipif(sys.platform == "win32", reason="/dev/fd is POSIX only")
def test_os_pipe_path(plain) -> None:  # type: ignore[no-untyped-def]
    r, w = os.pipe()
    feeder = threading.Thread(target=lambda: (os.write(w, TEXT.encode()), os.close(w)))
    feeder.start()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with PeffReader(f"/dev/fd/{r}") as reader:
                assert (reader.header, list(reader)) == plain
    finally:
        feeder.join(timeout=10)
        os.close(r)


@pytest.mark.skipif(sys.platform == "win32", reason="/dev/fd is POSIX only")
@pytest.mark.parametrize("suffix", ["", *COMPRESSORS])
def test_pipe_short_first_chunk(plain, suffix: str) -> None:  # type: ignore[no-untyped-def]
    # The first read of a pipe returns only what the writer has sent so far; the
    # format sniff must wait for the whole magic number, not decide on one byte.
    data = COMPRESSORS[suffix](TEXT.encode()) if suffix else TEXT.encode()
    r, w = os.pipe()

    def feed() -> None:
        os.write(w, data[:1])
        time.sleep(0.2)
        os.write(w, data[1:])
        os.close(w)

    feeder = threading.Thread(target=feed)
    feeder.start()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert read_peff(f"/dev/fd/{r}") == plain
    finally:
        feeder.join(timeout=10)
        os.close(r)


_NO_LZMA_BZ2 = """
import sys, warnings
sys.modules["_lzma"] = None
sys.modules["_bz2"] = None
import pefftacular
from pefftacular import PeffError, read_peff
warnings.simplefilter("ignore")
header, entries = read_peff(sys.argv[1])
assert entries
for path, module in ((sys.argv[2], "lzma"), (sys.argv[3], "bz2")):
    try:
        read_peff(path)
    except PeffError as e:
        assert module in str(e), e
    else:
        raise SystemExit(f"{path} read without {module}")
print("ok")
"""


def test_import_without_lzma_and_bz2(tmp_path: Path) -> None:
    # Minimal Python builds (pyenv, slim images) can lack _lzma and _bz2.
    plain_path = tmp_path / "p.peff"
    plain_path.write_bytes(TEXT.encode())
    xz = tmp_path / "x.peff.xz"
    xz.write_bytes(lzma.compress(TEXT.encode()))
    bz = tmp_path / "x.peff.bz2"
    bz.write_bytes(bz2.compress(TEXT.encode()))
    result = subprocess.run(
        [sys.executable, "-c", _NO_LZMA_BZ2, str(plain_path), str(xz), str(bz)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"


@pytest.mark.parametrize("compress", [False, True])
def test_non_utf8_input_raises_peff_parse_error(tmp_path: Path, compress: bool) -> None:
    data = "# PEFF 1.0\n# //\n# Prefix=sp\n# //\n>sp:P1 \\PName=caf\xe9\nMK\n".encode("latin-1")
    path = tmp_path / ("latin1.peff.gz" if compress else "latin1.peff")
    path.write_bytes(gzip.compress(data) if compress else data)
    with pytest.raises(PeffParseError, match="Cannot read the input") as info:
        read_peff(path)
    assert isinstance(info.value.__cause__, UnicodeDecodeError)
