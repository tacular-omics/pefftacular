"""SequenceEntry.from_fasta / to_fasta / to_proforma."""

from __future__ import annotations

import dataclasses
import re
import string
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from pefftacular import (
    ModRes,
    ModResPsi,
    ModResUnimod,
    PeffError,
    SequenceEntry,
    VariantSimple,
    read_peff,
)

FIXTURES = Path(__file__).parent / "fixtures"
UNIPROT = "sp|P31946|1433B_HUMAN 14-3-3 protein beta/alpha OS=Homo sapiens OX=9606 GN=YWHAB PE=1 SV=3"

# ---------------------------------------------------------------------------
# from_fasta / to_fasta
# ---------------------------------------------------------------------------


def test_from_fasta_uniprot() -> None:
    e = SequenceEntry.from_fasta(">" + UNIPROT, "MTMDK\nSELVQ K\n")
    assert (e.prefix, e.db_unique_id, e.id) == ("sp", "P31946", "1433B_HUMAN")
    assert e.pname == "14-3-3 protein beta/alpha"
    assert (e.tax_name, e.ncbi_tax_id, e.gname, e.pe, e.sv) == ("Homo sapiens", 9606, "YWHAB", 1, 3)
    assert e.sequence == "MTMDKSELVQK" and e.length == 11
    assert e.extra == {}
    assert e.to_fasta() == (UNIPROT, "MTMDKSELVQK")


@pytest.mark.parametrize(
    ("header", "kwargs", "expected"),
    [
        ("gi|12345|ref|NP_000001.1| some protein", {}, ("gi", "12345", None, "some protein")),
        ("nxp:NX_P01308-1 Insulin", {}, ("nxp", "NX_P01308-1", None, "Insulin")),
        ("ENSP0001\tthing X=1", {"prefix": "ens"}, ("ens", "ENSP0001", None, "thing")),
        ("sp|P1|N_HUMAN", {"prefix": "up"}, ("up", "P1", "N_HUMAN", None)),
        ("nr:gi|136429|sp|P00761.1|TRYP_PIG x", {}, ("nr", "gi|136429|sp|P00761.1|TRYP_PIG", None, "x")),
        ("sp|A:B|N_HUMAN", {}, ("sp", "A:B", "N_HUMAN", None)),
    ],
)
def test_from_fasta_identifier_shapes(header: str, kwargs: dict, expected: tuple) -> None:
    e = SequenceEntry.from_fasta(header, "MK", **kwargs)
    assert (e.prefix, e.db_unique_id, e.id, e.pname) == expected


def test_from_fasta_extra_and_bad_ints() -> None:
    e = SequenceEntry.from_fasta("sp|P1|N_HUMAN Name OX=abc PE=1 PE=2 XY=z w", "MK")
    assert e.ncbi_tax_id is None and e.pe == 1
    assert e.extra == {"OX": "abc", "PE": "2", "XY": "z w"}


@pytest.mark.parametrize("header", ["", ">", "   "])
def test_from_fasta_empty_header(header: str) -> None:
    with pytest.raises(PeffError, match="Empty"):
        SequenceEntry.from_fasta(header, "MK")


def test_from_fasta_needs_prefix() -> None:
    with pytest.raises(PeffError, match="prefix") as info:
        SequenceEntry.from_fasta("ENSP0001 thing", "MK")
    assert any("prefix=" in n for n in info.value.__notes__)


def test_to_fasta_drops_annotations_and_writes_minimal_header() -> None:
    e = SequenceEntry(prefix="nxp", db_unique_id="NX_1", sequence="MKV", mod_res_psi=(ModResPsi((1,), "MOD:1", "x"),))
    assert e.to_fasta() == ("nxp|NX_1", "MKV")


@pytest.mark.parametrize("path", sorted(FIXTURES.glob("*_Valid.peff")) + sorted(FIXTURES.glob("[a-z]*.peff")))
def test_fixture_entries_round_trip_their_fasta_fields(path: Path) -> None:
    _, entries = read_peff(path)
    for entry in entries:
        back = SequenceEntry.from_fasta(*entry.to_fasta())
        if "|" in entry.db_unique_id:
            assert entry.to_fasta()[0].startswith(f"{entry.prefix}:{entry.db_unique_id}")
            assert back.id is None
            entry = dataclasses.replace(entry, id=None)
        for name in ("prefix", "db_unique_id", "id", "tax_name", "ncbi_tax_id", "gname", "pe", "sv", "sequence"):
            assert getattr(back, name) == getattr(entry, name)
        if entry.pname and not re.search(r"\s[A-Za-z_]\w*=", " " + entry.pname):
            assert back.pname == entry.pname


def test_to_fasta_header_reads_back_in_fastatacular() -> None:
    fastatacular = pytest.importorskip("fastatacular")
    import io

    e = SequenceEntry.from_fasta(UNIPROT, "MTMDK")
    header, seq = e.to_fasta()
    [fe] = fastatacular.read_fasta(io.StringIO(f">{header}\n{seq}\n"))
    assert (fe.prefix, fe.accession, fe.entry_name, fe.pname) == ("sp", "P31946", "1433B_HUMAN", e.pname)
    assert (fe.os_name, fe.ncbi_tax_id, fe.gname, fe.pe, fe.sv) == ("Homo sapiens", 9606, "YWHAB", 1, 3)
    assert SequenceEntry.from_fasta(fe.raw_header, fe.sequence) == e


_word = st.text(string.ascii_letters + string.digits + "-_.,()/'", min_size=1, max_size=8)
_text = st.lists(_word, min_size=1, max_size=4).map(" ".join)
_token = st.text(string.ascii_letters + string.digits + "_.-", min_size=1, max_size=10)


@given(
    prefix=_token,
    acc=_token,
    name=st.none() | _token,
    pname=st.none() | _text,
    tax=st.none() | _text,
    ox=st.none() | st.integers(0, 10**7),
    gn=st.none() | _word,
    pe=st.none() | st.integers(1, 5),
    sv=st.none() | st.integers(0, 99),
    seq=st.text("ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=40),
)
def test_fasta_round_trip_property(prefix, acc, name, pname, tax, ox, gn, pe, sv, seq) -> None:  # type: ignore[no-untyped-def]
    e = SequenceEntry(
        prefix=prefix,
        db_unique_id=acc,
        id=name,
        sequence=seq,
        length=len(seq),
        pname=pname,
        tax_name=tax,
        ncbi_tax_id=ox,
        gname=gn,
        pe=pe,
        sv=sv,
    )
    assert SequenceEntry.from_fasta(*e.to_fasta()) == e


# ---------------------------------------------------------------------------
# to_proforma
# ---------------------------------------------------------------------------


def _parse_proforma(s: str) -> tuple[str, dict[int, list[str]], dict[str, int]]:
    """Tiny inverse of to_proforma for the subset it writes."""
    unknown: dict[str, int] = {}
    if "?" in s:
        head, s = s.split("?", 1)
        for m in re.finditer(r"\[([^\]]+)\](?:\^(\d+))?", head):
            unknown[m[1]] = int(m[2] or 1)
    seq = ""
    at: dict[int, list[str]] = {}
    i = 0
    while i < len(s):
        if s[i] == "[":
            j = s.index("]", i)
            at.setdefault(len(seq), []).append(s[i + 1 : j])
            i = j + 1
        else:
            seq += s[i]
            i += 1
    return seq, at, unknown


def test_to_proforma_basic() -> None:
    e = SequenceEntry(
        prefix="sp",
        db_unique_id="P1",
        sequence="MSTKSY",
        mod_res_psi=(ModResPsi((2, 5), "MOD:00046", "O-phospho-L-serine"), ModResPsi(("?", "?"), "MOD:00047", "x")),
        mod_res_unimod=(ModResUnimod((1,), "UNIMOD:1", "Acetyl"), ModResUnimod((2,), "21", "Phospho")),
        mod_res=(
            ModRes((6,), "MOD:00048", "tyr"),
            ModRes((4,), "", "N-linked (GlcNAc...)"),
            ModRes((3,), "UNIMOD:21", "Phospho"),
        ),
    )
    assert e.to_proforma() == "[MOD:00047]^2?MS[MOD:00046]TKS[MOD:00046]Y[MOD:00048]"
    assert e.to_proforma(mods="unimod") == "M[UNIMOD:1]S[UNIMOD:21]T[UNIMOD:21]KSY"


def test_to_proforma_names_when_no_accession() -> None:
    e = SequenceEntry(
        prefix="x", db_unique_id="1", sequence="MK", mod_res_unimod=(ModResUnimod((2, "?"), "", "Methyl"),)
    )
    assert e.to_proforma(mods="unimod") == "[U:Methyl]?MK[U:Methyl]"
    e = dataclasses.replace(e, mod_res_psi=(ModResPsi((1,), "", "methylated"),))
    assert e.to_proforma() == "M[M:methylated]K"


def test_to_proforma_variants() -> None:
    e = SequenceEntry(
        prefix="sp",
        db_unique_id="P1",
        sequence="MSTKSY",
        variant_simple=(VariantSimple(2, "A"), VariantSimple(5, "*"), VariantSimple(3, "S")),
        mod_res_psi=(ModResPsi((2, 3, 6), "MOD:00046", "p"),),
    )
    assert e.to_proforma(variants=e.variant_simple[:1]) == "MAT[MOD:00046]KSY[MOD:00046]"
    assert e.to_proforma(variants=e.variant_simple[1:2]) == "MS[MOD:00046]T[MOD:00046]K"
    assert e.to_proforma(variants=[e.variant_simple[2]]) == "MS[MOD:00046]SKSY[MOD:00046]"
    assert e.to_proforma(variants=[VariantSimple(1, "*")]) == ""


@pytest.mark.parametrize(
    ("kwargs", "entry_kwargs", "match"),
    [
        ({"mods": "resid"}, {}, "mods must be"),
        ({"variants": [VariantSimple(9, "A")]}, {}, "VariantSimple position 9"),
        ({"variants": [VariantSimple("x", "A")]}, {}, "VariantSimple position 'x'"),
        ({"variants": [VariantSimple(1, "AA")]}, {}, "single letter"),
        ({"variants": [VariantSimple(1, "A"), VariantSimple(1, "C")]}, {}, "Two different"),
        ({}, {"mod_res_psi": (ModResPsi((0,), "MOD:1", "x"),)}, "ModResPsi position 0"),
        ({}, {"mod_res_psi": (ModResPsi((1,), "MOD:[1]", "x"),)}, "square bracket"),
    ],
)
def test_to_proforma_errors(kwargs: dict, entry_kwargs: dict, match: str) -> None:
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence="MKV", **entry_kwargs)
    with pytest.raises(PeffError, match=match):
        e.to_proforma(**kwargs)


def test_to_proforma_same_variant_twice_is_fine() -> None:
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence="MKV")
    assert e.to_proforma(variants=[VariantSimple(2, "R"), VariantSimple(2, "R")]) == "MRV"


@given(
    seq=st.text("ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=30),
    data=st.data(),
    mods=st.sampled_from(["psimod", "unimod"]),
)
def test_proforma_round_trip_property(seq: str, data: st.DataObject, mods: str) -> None:
    n = len(seq)
    pos = st.integers(1, n) | st.just("?")
    acc = st.integers(1, 99999).map(str)
    psi = data.draw(st.lists(st.tuples(st.lists(pos, min_size=1, max_size=3), acc), max_size=4))
    uni = data.draw(st.lists(st.tuples(st.lists(pos, min_size=1, max_size=3), acc), max_size=4))
    e = SequenceEntry(
        prefix="x",
        db_unique_id="1",
        sequence=seq,
        mod_res_psi=tuple(ModResPsi(tuple(p), f"MOD:{a}", "n") for p, a in psi),
        mod_res_unimod=tuple(ModResUnimod(tuple(p), f"UNIMOD:{a}", "n") for p, a in uni),
    )
    source, cv = (psi, "MOD") if mods == "psimod" else (uni, "UNIMOD")
    expected_at: dict[int, list[str]] = {}
    expected_unknown: dict[str, int] = {}
    for positions, a in source:
        for p in positions:
            if p == "?":
                expected_unknown[f"{cv}:{a}"] = expected_unknown.get(f"{cv}:{a}", 0) + 1
            else:
                expected_at.setdefault(p, []).append(f"{cv}:{a}")
    got_seq, got_at, got_unknown = _parse_proforma(e.to_proforma(mods=mods))  # type: ignore[arg-type]
    assert got_seq == seq
    assert got_at == expected_at
    assert got_unknown == expected_unknown


@pytest.mark.parametrize("path", sorted(FIXTURES.glob("*.peff")))
def test_every_valid_fixture_entry_renders(path: Path) -> None:
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            _, entries = read_peff(path)
        except PeffError:
            pytest.skip("invalid fixture")
    for e in entries:
        for mods in ("psimod", "unimod"):
            try:
                s = e.to_proforma(mods=mods)  # type: ignore[arg-type]
            except PeffError:
                continue  # fixtures with out-of-range positions (INValid files)
            assert _parse_proforma(s)[0] == e.sequence
