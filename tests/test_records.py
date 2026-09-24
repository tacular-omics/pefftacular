import io
import re
import warnings
from pathlib import Path

import pytest

from pefftacular import RECORD_KEYS, PeffReader, SequenceEntry, read_peff, to_records

FIXTURES = Path(__file__).parent / "fixtures"
README = Path(__file__).resolve().parents[1] / "README.md"
PARSEABLE = [
    "PEFF_AnnotID_Insulin_Valid.peff",
    "PEFF_Minimal_Valid.peff",
    "PEFF_Tiny_Valid.peff",
    "SmallTestDB-PEFF1.0.peff",
    "UniProtExport_3prot.peff",
    "complex.peff",
    "custom_keys.peff",
    "minimal.peff",
    "multidb.peff",
    "proteoform_ENST00000000412.peff",
]

# record key -> PEFF description key, in the writer's order
_PEFF_KEYS = {
    "id": "ID",
    "db_unique_id_key": "DbUniqueId",
    "length": "Length",
    "pname": "PName",
    "gname": "GName",
    "ncbi_tax_id": "NcbiTaxId",
    "tax_name": "TaxName",
    "sv": "SV",
    "ev": "EV",
    "pe": "PE",
    "decoy": "Decoy",
    "comment": "Comment",
    "variant_simple": "VariantSimple",
    "variant_complex": "VariantComplex",
    "mod_res_unimod": "ModResUnimod",
    "mod_res_psi": "ModResPsi",
    "mod_res": "ModRes",
    "processed": "Processed",
    "disulfide_bond": "DisulfideBond",
    "proteoform": "Proteoform",
}


def _escape(s: str) -> str:
    # the writer's free-text escaping (see _writer._escape_component) for the fixtures' values
    from pefftacular._writer import _escape_component

    return _escape_component(s)


def _description_line(rec: dict) -> str:
    parts = [f">{rec['prefix']}:{rec['db_unique_id']}"]
    for key, peff in _PEFF_KEYS.items():
        value = rec[key]
        if value is None:
            continue
        if isinstance(value, bool):
            value = "true" if value else "false"
        elif key in ("pname", "gname", "tax_name", "comment"):
            value = _escape(value)
        parts.append(f"\\{peff}={value}")
    for key in ("custom_values", "extra"):
        if rec[key] is not None:
            parts.append(rec[key])
    return " ".join(parts)


def _read(path: Path | str | io.StringIO):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_peff(path)


def _records(source):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return to_records(source)


@pytest.mark.parametrize("name", PARSEABLE)
def test_keys_and_roundtrip(name: str) -> None:
    path = FIXTURES / name
    _, entries = _read(path)
    records = _records(path)
    assert len(records) == len(entries)
    header_text = "".join(
        line for line in path.read_text(encoding="utf-8-sig").splitlines(True) if line.startswith("#")
    )
    for rec, entry in zip(records, entries, strict=True):
        assert tuple(rec) == RECORD_KEYS
        for key in RECORD_KEYS:
            value = getattr(entry, key)
            if key in ("custom_values", "extra") or isinstance(value, tuple):
                assert (rec[key] is None) == (not value)
            else:
                assert rec[key] == value
        # the flat record carries everything: rebuild the entry from it and re-parse
        text = header_text + _description_line(rec) + "\n" + rec["sequence"] + "\n"
        _, (again,) = _read(io.StringIO(text))
        assert again == entry


def test_sources_agree() -> None:
    path = FIXTURES / "SmallTestDB-PEFF1.0.peff"
    expected = _records(path)
    assert _records(str(path)) == expected
    assert _records(io.StringIO(path.read_text(encoding="utf-8"))) == expected
    entries = _read(path)[1]
    assert _records(entries) == expected
    assert [e.to_record() for e in entries] == expected
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with PeffReader(path) as reader:
            assert reader.to_records() == expected


def test_custom_values_use_header_defs() -> None:
    rec = _records(FIXTURES / "custom_keys.peff")
    assert any(r["custom_values"] and r["custom_values"].startswith("\\") for r in rec)


def test_values() -> None:
    rec = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="MK", decoy=True, extra={"Foo": "bar"}).to_record()
    assert rec["decoy"] is True
    assert rec["extra"] == "\\Foo=bar"
    assert rec["variant_simple"] is None and rec["length"] is None


def test_empty() -> None:
    assert _records(io.StringIO("# PEFF 1.0\n# //\n")) == []
    assert to_records([]) == []


def test_empty_file(tmp_path: Path) -> None:
    p = tmp_path / "empty.peff"
    p.write_text("# PEFF 1.0\n")
    assert _records(p) == []


def test_pandas() -> None:
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(_records(FIXTURES / "SmallTestDB-PEFF1.0.peff"))
    assert list(df.columns) == list(RECORD_KEYS)
    assert df["variant_simple"].notna().sum() == 9


def test_polars() -> None:
    pl = pytest.importorskip("polars")
    df = pl.DataFrame(_records(FIXTURES / "SmallTestDB-PEFF1.0.peff"))
    assert df.columns == list(RECORD_KEYS)
    assert df.filter(pl.col("variant_simple").is_not_null()).height == 9


def _readme_block() -> str:
    section = README.read_text(encoding="utf-8").split("## Tables with pandas or polars", 1)[1]
    return re.search(r"```python\n(.*?)```", section, re.S).group(1)  # type: ignore[union-attr]


def test_readme_example(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("pandas")
    pytest.importorskip("polars")
    (tmp_path / "proteins.peff").write_bytes((FIXTURES / "SmallTestDB-PEFF1.0.peff").read_bytes())
    monkeypatch.chdir(tmp_path)
    ns: dict[str, object] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(_readme_block(), ns)
    assert len(ns["with_variants"]) == 9  # type: ignore[arg-type]
