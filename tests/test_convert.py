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
        ("gi|12345|ref|NP_000001.1| some protein", {}, ("ref", "NP_000001.1", None, "some protein")),
        ("gi|4557757|ref|NP_000240.1| desc", {}, ("ref", "NP_000240.1", None, "desc")),
        ("gi|7|gb|AAA12345.1|", {}, ("gb", "AAA12345.1", None, None)),
        ("gi|7|emb|CAA12345.1|", {}, ("emb", "CAA12345.1", None, None)),
        ("gi|7|dbj|BAA12345.1|", {}, ("dbj", "BAA12345.1", None, None)),
        ("gi|7|pdb|1ABC|A", {}, ("pdb", "1ABC_A", None, None)),
        ("gi|1|pdb|1MBA|B", {}, ("pdb", "1MBA_B", None, None)),
        ("gi|7|pdb|1ABC|", {}, ("pdb", "1ABC", None, None)),
        ("DECOY_gi|1|ref|NP_2.1|", {}, ("DECOY_ref", "NP_2.1", None, None)),
        ("rev_gi|1|ref|NP_2.1| d", {}, ("rev_ref", "NP_2.1", None, "d")),
        ("CONTAM_gi|1|sp|P1.1|N_HUMAN", {}, ("CONTAM_sp", "P1.1", "N_HUMAN", None)),
        ("gi|1|ref|NP_3.1|x|y", {}, ("ref", "NP_3.1", None, None)),
        ("gi|1|sp|P1.1|N_HUMAN|extra", {}, ("sp", "P1.1", "N_HUMAN", None)),
        ("fungi|1|ref|NP_2.1|", {}, ("fungi", "1", None, None)),
        ("fungi|5|ab|Q|", {}, ("fungi", "5", None, None)),
        ("decoy_gi|1|ref|NP_2.1|", {}, ("decoy_ref", "NP_2.1", None, None)),
        ("REV-2-gi|1|ref|NP_2.1|", {}, ("REV-2-ref", "NP_2.1", None, None)),
        ("XXX_gi|1|ref|NP_2.1|", {}, ("XXX_gi", "1", None, None)),
        ("gi|136429|sp|P00761.1|TRYP_PIG x", {}, ("sp", "P00761.1", "TRYP_PIG", "x")),
        ("gi|7|tr|Q12345|", {}, ("tr", "Q12345", None, None)),
        ("gi|7|ref|NP_1.1| p", {"prefix": "ncbi"}, ("ncbi", "NP_1.1", None, "p")),
        ("gi|12345 lone gi", {}, ("gi", "12345", None, "lone gi")),
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


# Valid fixtures only. Filter by name in Python: glob is case-insensitive on Windows.
_VALID_FIXTURES = [p for p in sorted(FIXTURES.glob("*.peff")) if p.name.endswith("_Valid.peff") or p.name[0].islower()]


@pytest.mark.parametrize("path", _VALID_FIXTURES, ids=lambda p: p.name)
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
        ({"variants": [VariantSimple(1, "A"), VariantSimple(1, "C")]}, {}, "two different"),
        ({}, {"mod_res_psi": (ModResPsi((0,), "MOD:1", "x"),)}, "ModResPsi position 0"),
        ({}, {"mod_res_psi": (ModResPsi((1,), "MOD:[1]", "x"),)}, "square bracket"),
    ],
)
def test_to_proforma_errors(kwargs: dict, entry_kwargs: dict, match: str) -> None:
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence="MKV", **entry_kwargs)
    with pytest.raises(PeffError, match=match) as info:
        e.to_proforma(**kwargs)
    if "mods" not in kwargs:
        assert str(info.value).startswith("x:1: ")  # names the entry
        assert e.to_proforma(**kwargs, errors="skip") is None
    else:
        with pytest.raises(PeffError, match=match):  # a bad argument raises even with skip
            e.to_proforma(**kwargs, errors="skip")


def test_to_proforma_site_past_the_end_names_accession_and_position() -> None:
    # As in the neXtProt human PEFF, where 12 entries list sites past the sequence end.
    e = SequenceEntry(
        prefix="nxp", db_unique_id="NX_P1-2", sequence="MKV", mod_res_psi=(ModResPsi((2, 7), "MOD:00046", "p"),)
    )
    with pytest.raises(PeffError, match=r"^nxp:NX_P1-2: ModResPsi position 7 is not a residue position \(1\.\.3\)"):
        e.to_proforma()
    ok = SequenceEntry(prefix="nxp", db_unique_id="NX_P2", sequence="MKV")
    forms = [p for x in (e, ok) if (p := x.to_proforma(errors="skip")) is not None]
    assert forms == ["MKV"]
    assert ok.to_proforma(errors="skip") == ok.to_proforma() == "MKV"
    with pytest.raises(PeffError, match="errors must be"):
        ok.to_proforma(errors="ignore")  # type: ignore[call-overload]


def test_to_proforma_same_site_in_psi_and_generic_written_once() -> None:
    e = SequenceEntry(
        prefix="x",
        db_unique_id="1",
        sequence="MSTK",
        mod_res_psi=(ModResPsi((2, "?"), "MOD:00046", "p"),),
        mod_res=(ModRes((2, 3, "?"), "MOD:00046", "p"),),
    )
    assert e.to_proforma() == "[MOD:00046]?MS[MOD:00046]T[MOD:00046]K"


def test_to_proforma_unknown_sites_dropped_after_truncation() -> None:
    e = SequenceEntry(
        prefix="x", db_unique_id="1", sequence="MSTK", mod_res_psi=(ModResPsi((2, "?"), "MOD:00046", "p"),)
    )
    assert e.to_proforma() == "[MOD:00046]?MS[MOD:00046]TK"
    assert e.to_proforma(variants=[VariantSimple(3, "*")]) == "MS[MOD:00046]"


def test_dataclass_field_types_are_types_not_strings() -> None:
    # Introspection (cattrs, pydantic dataclasses) needs real annotation objects.
    import pefftacular._models as models

    for cls in vars(models).values():
        if dataclasses.is_dataclass(cls) and isinstance(cls, type):
            for f in dataclasses.fields(cls):
                assert not isinstance(f.type, str), (cls.__name__, f.name, f.type)


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
            elif f"{cv}:{a}" not in expected_at.setdefault(p, []):
                expected_at[p].append(f"{cv}:{a}")
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


@pytest.mark.parametrize("seq", ["MK*V", "MK-V", "MKvV", "MK1V"])
def test_to_proforma_residue_outside_a_to_z_raises(seq: str) -> None:
    # PEFF allows "*" (interruption) and "-" (gap); a ProForma residue is one of A-Z.
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence=seq)
    with pytest.raises(PeffError, match=r"^x:1: residue .* at position 3 cannot be written as ProForma"):
        e.to_proforma()
    assert e.to_proforma(errors="skip") is None
    assert e.to_proforma(variants=[VariantSimple(3, "*")]) == "MK"  # truncated before it: fine


def test_to_proforma_same_modification_in_two_spellings_is_written_once() -> None:
    # \ModResPsi and \ModRes both list MOD:00046 on K3, in different case and spacing.
    e = SequenceEntry(
        prefix="x",
        db_unique_id="1",
        sequence="PEK",
        mod_res_psi=(ModResPsi(positions=(3, "?"), accession="MOD:00046", name="O-phospho-L-serine"),),
        mod_res=(ModRes(positions=(3, "?"), accession=" mod:00046", name="O-phospho-L-serine"),),
    )
    assert e.to_proforma() == "[MOD:00046]?PEK[MOD:00046]"
    pt = pytest.importorskip("peptacular")
    site_only = dataclasses.replace(
        e, mod_res_psi=(ModResPsi(positions=(3,), accession="MOD:00046", name="p"),), mod_res=()
    )
    one = site_only.to_proforma()
    assert one == "PEK[MOD:00046]"
    two = dataclasses.replace(site_only, mod_res=(ModRes(positions=(3,), accession="mod:00046", name="p"),))
    assert pt.mass(two.to_proforma()) == pytest.approx(pt.mass(one))
    assert pt.mass(one) == pytest.approx(pt.mass("PEK") + 79.966331, abs=1e-4)


@pytest.mark.parametrize("new", ["k", "\u00e9", "1"])
def test_to_proforma_variant_outside_a_to_z_names_the_variant(new: str) -> None:
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence="MKV")
    with pytest.raises(PeffError, match=r"^x:1: VariantSimple at 2: new amino acid .* is not a single letter A-Z"):
        e.to_proforma(variants=[VariantSimple(2, new)])


@pytest.mark.parametrize(("accession", "name"), [("", "a|b"), ("", "a#b"), ("MOD:00046|INFO:x", "p"), ("MOD:1#g", "p")])
def test_to_proforma_pipe_or_hash_in_modification_raises(accession: str, name: str) -> None:
    # No ProForma escape: "|" would start a second tag and "#" a group label.
    e = SequenceEntry(prefix="x", db_unique_id="1", sequence="MKV", mod_res_psi=(ModResPsi((2,), accession, name),))
    with pytest.raises(PeffError, match=r"^x:1: modification .* contains '[|#]' and cannot be written as ProForma"):
        e.to_proforma()
    assert e.to_proforma(errors="skip") is None


# ---------------------------------------------------------------------------
# to_proforma output read by a real ProForma parser (peptacular)
# ---------------------------------------------------------------------------

# Mod names: no square bracket, "|" or "#" (to_proforma refuses them, tested above), no
# leading/trailing space (a ProForma reader may trim it) and a leading letter, as in every
# PSI-MOD and Unimod name ("[M:+1]" would be a mass shift). "?" and parens are kept on purpose.
_MOD_NAME = st.text(string.ascii_letters + string.digits + " -_.,:;()'?+*^", min_size=1, max_size=12).filter(
    lambda s: s == s.strip() and s[0].isalpha()
)
_CV = {"psimod": ("MOD", "M"), "unimod": ("UNIMOD", "U")}


def _accession(cv: str) -> st.SearchStrategy[str]:
    number = st.integers(1, 99999)
    return st.one_of(
        number.map(lambda n: f"{cv}:{n:05d}" if cv == "MOD" else f"{cv}:{n}"),
        number.map(str),  # bare number: to_proforma adds the CV prefix
        st.just(""),  # no accession: written by name
    )


def _tag_key(tag: object) -> tuple[str, str, str]:
    """(kind, CV, value) of a peptacular tag, or of a ``"CV:value"`` string (CV case-insensitive)."""
    if isinstance(tag, str):
        cv, _, value = tag.partition(":")
        return ("name", {"M": "MOD", "U": "UNIMOD"}[cv], value) if cv in ("M", "U") else ("acc", cv.upper(), value)
    kind = "acc" if hasattr(tag, "accession") else "name"
    return kind, tag.cv.value, tag.accession if kind == "acc" else tag.name  # type: ignore[attr-defined]


def _counts(mods: object) -> dict[tuple[str, str, str], int]:
    out: dict[tuple[str, str, str], int] = {}
    for m in mods.mods:  # type: ignore[attr-defined]
        (tag,) = m.value.tags  # one tag per modification: a second would mean a lost "|"
        key = _tag_key(tag)
        assert key not in out, f"{key} written twice"  # e.g. MOD:00046 and mod:00046: mass counted twice
        out[key] = m.count
    return out


@st.composite
def _proforma_case(draw: st.DrawFn) -> tuple[SequenceEntry, str, tuple[VariantSimple, ...]]:
    vocab = draw(st.sampled_from(sorted(_CV)))
    # Mostly A-Z; "*", "-" and lowercase letters cannot be written and must raise when kept.
    residue = st.sampled_from("ACDEFGHIKLMNPQRSTVWYBJOUXZ") | st.sampled_from("*-acdkmz")
    seq = draw(st.text(residue, min_size=1, max_size=25))
    n = len(seq)
    # Termini on purpose as well as anywhere, and unknown sites.
    pos = st.sampled_from([1, n]) | st.integers(1, n) | st.just("?")
    positions = st.lists(pos, min_size=1, max_size=3).map(tuple)

    def mods(cls: type, cv: str) -> st.SearchStrategy:
        return st.lists(st.builds(cls, positions=positions, accession=_accession(cv), name=_MOD_NAME), max_size=3)

    variants = draw(
        st.lists(
            st.builds(VariantSimple, position=st.integers(1, n), new_amino_acid=st.sampled_from([*"ACDKRSTY*"])),
            max_size=3,
            unique_by=lambda v: v.position,
        )
    )
    entry = SequenceEntry(
        prefix="x",
        db_unique_id="1",
        sequence=seq,
        mod_res_psi=tuple(draw(mods(ModResPsi, "MOD"))),
        mod_res_unimod=tuple(draw(mods(ModResUnimod, "UNIMOD"))),
        # \ModRes from either vocabulary or another one: only the matching CV is written.
        mod_res=tuple(
            draw(
                st.lists(
                    st.builds(
                        ModRes,
                        positions=positions,
                        accession=st.builds(
                            lambda acc, case, before, after: before + case(acc) + after,
                            _accession("MOD").filter(bool)
                            | _accession("UNIMOD").filter(lambda a: ":" in a)
                            | st.just("RESID:AA0037"),
                            # The CV prefix is matched after .strip().upper().
                            st.sampled_from([str, str.lower, str.title]),
                            st.sampled_from(["", " "]),
                            st.sampled_from(["", " "]),
                        ),
                        name=_MOD_NAME,
                    ),
                    max_size=3,
                )
            )
        ),
        variant_simple=tuple(variants),
    )
    applied = tuple(draw(st.lists(st.sampled_from(variants), unique=True))) if variants else ()
    return entry, vocab, applied


def _expected(
    entry: SequenceEntry, vocab: str, variants: tuple[VariantSimple, ...]
) -> tuple[str, dict[int, dict], dict]:
    """The residues, per-position tags and unknown-site counts to_proforma should write.

    Tags are keyed by their normalised (kind, CV, value): "MOD:1" and " mod:1" are one
    modification and must be written once per site.
    """
    cv, short = _CV[vocab]
    seq = list(entry.sequence)
    stops = [v.position for v in variants if v.new_amino_acid == "*"]
    end = min(stops, default=len(seq) + 1) - 1
    subs = {v.position: v.new_amino_acid for v in variants if v.new_amino_acid != "*"}
    for p, aa in subs.items():
        seq[p - 1] = aa

    def tag(m: ModRes | ModResPsi | ModResUnimod) -> tuple[str, str, str]:
        acc = m.accession.strip()
        if not acc:
            return ("name", cv, m.name)
        return _tag_key(acc if ":" in acc else f"{cv}:{acc}")

    own = entry.mod_res_psi if vocab == "psimod" else entry.mod_res_unimod
    generic = [m for m in entry.mod_res if m.accession.strip().upper().startswith(cv + ":")]
    at: dict[int, dict] = {}
    unknown: dict = {}
    for source in (own, generic):
        source_unknown: dict = {}
        for m in source:
            for p in m.positions:
                if p == "?":
                    source_unknown[tag(m)] = source_unknown.get(tag(m), 0) + 1
                elif p not in subs and p <= end:
                    at.setdefault(p - 1, {})[tag(m)] = 1
        for t, c in source_unknown.items():
            unknown[t] = max(unknown.get(t, 0), c)
    if end < len(seq):
        unknown = {}
    return "".join(seq[:end]), at, unknown


@given(case=_proforma_case())
def test_to_proforma_reads_back_in_peptacular(case: tuple[SequenceEntry, str, tuple[VariantSimple, ...]]) -> None:
    pt = pytest.importorskip("peptacular")
    entry, vocab, variants = case
    expected_seq, expected_at, expected_unknown = _expected(entry, vocab, variants)

    # Raises if and only if a residue outside A-Z survives substitution and truncation.
    if any(not ("A" <= aa <= "Z") for aa in expected_seq):
        with pytest.raises(PeffError, match="cannot be written as ProForma"):
            entry.to_proforma(mods=vocab, variants=variants)  # type: ignore[arg-type]
        assert entry.to_proforma(mods=vocab, variants=variants, errors="skip") is None  # type: ignore[arg-type]
        return
    proforma = entry.to_proforma(mods=vocab, variants=variants)  # type: ignore[arg-type]
    annotation = pt.parse(proforma)
    assert annotation.stripped_sequence == expected_seq
    assert {i: _counts(m) for i, m in annotation.internal_mods.items()} == expected_at
    assert _counts(annotation.unknown_mods) == expected_unknown
    # PEFF cannot mark a terminal modification: everything sits on a residue.
    assert not annotation.has_nterm_mods and not annotation.has_cterm_mods
