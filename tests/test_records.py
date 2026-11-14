import dataclasses
import io
import re
import warnings
from pathlib import Path

import pytest

from pefftacular import (
    RECORD_KEYS,
    PeffParseError,
    PeffReader,
    PeffWriteError,
    SequenceEntry,
    from_records,
    read_peff,
    to_records,
)

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


# ---------------------------------------------------------------------------
# from_records: records back to entries
# ---------------------------------------------------------------------------

# The PSI "INValid" tiny example and the pre-1.0 draft are not readable PEFF 1.0 at all,
# so there is nothing to round-trip; test_excluded_fixtures_are_unreadable proves it.
_UNREADABLE = {"PEFF_Tiny_INValid1.peff", "SmallTestDB-PEFF0.9.peff"}
_ROUND_TRIP = sorted(p for p in FIXTURES.rglob("*.peff") if p.name not in _UNREADABLE)


def _entry_record(**overrides: object) -> dict:
    """One record for a minimal entry, with *overrides* applied."""
    entry = SequenceEntry(prefix="sp", db_unique_id="P1", sequence="ACDEFGHIKL", pname="Alpha")
    record = to_records([entry])[0]
    record.update(overrides)  # type: ignore[arg-type]
    return record


class TestFromRecordsRoundTrip:
    @pytest.mark.parametrize("path", _ROUND_TRIP, ids=lambda p: p.name)
    def test_fixture_round_trips(self, path: Path) -> None:
        header, entries = _read(path)
        assert from_records(to_records(entries), header) == entries

    def test_excluded_fixtures_are_unreadable(self) -> None:
        # The exclusion list above must stay honest: each name is unreadable, not just awkward.
        for name in _UNREADABLE:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with pytest.raises(PeffParseError):
                    read_peff(FIXTURES / name)

    def test_the_sweep_actually_covers_custom_values(self) -> None:
        # Without this the headline case could drop out of the sweep unnoticed.
        carriers = {
            path.name: sum(len(e.custom_values) for e in _read(path)[1])
            for path in _ROUND_TRIP
            if any(e.custom_values for e in _read(path)[1])
        }
        assert carriers == {"custom_keys.peff": 2}

    def test_custom_values_survive_both_directions(self) -> None:
        header, entries = _read(FIXTURES / "custom_keys.peff")
        records = to_records(entries)
        assert records[0]["custom_values"] is not None
        back = from_records(records, header)
        assert [e.custom_values for e in back] == [e.custom_values for e in entries]
        assert all(not e.extra for e in back)  # nothing fell through to the unknown-key bucket

    def test_entries_can_be_written_again(self) -> None:
        from pefftacular import write_peff

        header, entries = _read(FIXTURES / "custom_keys.peff")
        buf = io.StringIO()
        write_peff(header, from_records(to_records(entries), header), buf)
        buf.seek(0)
        assert _read(buf)[1] == entries


class TestFromRecordsRejects:
    def test_missing_key(self) -> None:
        record = _entry_record()
        del record["custom_values"]
        with pytest.raises(PeffWriteError, match=r"missing key\(s\) 'custom_values'") as err:
            from_records([record])
        assert err.value.index == 0

    def test_unknown_key(self) -> None:
        with pytest.raises(PeffWriteError, match=r"unknown key\(s\) 'oops'"):
            from_records([_entry_record(oops="x")])

    def test_custom_values_without_a_header(self) -> None:
        records = to_records(_read(FIXTURES / "custom_keys.peff")[1])
        with pytest.raises(PeffWriteError, match=r"no CustomKeyDef declares"):
            from_records(records)

    def test_custom_values_with_a_header_that_lacks_the_key(self) -> None:
        header, entries = _read(FIXTURES / "custom_keys.peff")
        stripped = tuple(dataclasses.replace(db, custom_key_defs=()) for db in header.databases)
        bare = dataclasses.replace(header, databases=stripped)
        with pytest.raises(PeffWriteError, match=r"no CustomKeyDef declares"):
            from_records(to_records(entries), bare)

    def test_a_key_set_twice(self) -> None:
        with pytest.raises(PeffWriteError, match=r"extra repeats \\PName="):
            from_records([_entry_record(extra=r"\PName=other")])

    @pytest.mark.parametrize("value", [9606.5, "nine", True, [9606]])
    def test_a_bad_int(self, value: object) -> None:
        with pytest.raises(PeffWriteError, match=r"ncbi_tax_id="):
            from_records([_entry_record(ncbi_tax_id=value)])

    @pytest.mark.parametrize("value", ["yes", 2, 1.5])
    def test_a_bad_bool(self, value: object) -> None:
        with pytest.raises(PeffWriteError, match=r"decoy="):
            from_records([_entry_record(decoy=value)])

    @pytest.mark.parametrize("key", ["prefix", "db_unique_id", "sequence"])
    def test_a_missing_identity_value(self, key: str) -> None:
        with pytest.raises(PeffWriteError, match=rf"{key} is missing"):
            from_records([_entry_record(**{key: None})])

    @pytest.mark.parametrize("key", ["prefix", "db_unique_id", "sequence", "variant_simple", "extra"])
    def test_an_empty_text_value(self, key: str) -> None:
        with pytest.raises(PeffWriteError, match=rf"{key} (is empty|holds no)"):
            from_records([_entry_record(**{key: ""})])

    def test_text_that_is_not_what_the_writer_emits(self) -> None:
        # PEFF tolerates a space between items, the writer never writes one, and a silent
        # rewrite is the thing this function exists to prevent.
        with pytest.raises(PeffWriteError, match=r"reads back as '\(3\|A\)\(10\|G\)'"):
            from_records([_entry_record(variant_simple="(3|A) (10|G)")])

    def test_not_a_mapping(self) -> None:
        with pytest.raises(PeffWriteError, match=r"is not a mapping"):
            from_records([["prefix", "sp"]])  # type: ignore[list-item]


class TestFromRecordsValues:
    def test_an_empty_string_is_a_value_not_an_absence(self) -> None:
        kept = from_records([_entry_record(pname="")])[0]
        gone = from_records([_entry_record(pname=None)])[0]
        assert (kept.pname, gone.pname) == ("", None)

    @pytest.mark.parametrize("missing", [None, float("nan")])
    def test_frame_missing_markers_are_absences(self, missing: object) -> None:
        assert from_records([_entry_record(pname=missing)])[0].pname is None

    @pytest.mark.parametrize("value", [9606, 9606.0, "9606", " 9606 "])
    def test_int_columns_accept_frame_values(self, value: object) -> None:
        assert from_records([_entry_record(ncbi_tax_id=value)])[0].ncbi_tax_id == 9606

    @pytest.mark.parametrize(
        ("value", "expected"),
        [(True, True), (False, False), (1, True), (0, False), ("true", True), ("false", False)],
    )
    def test_decoy_accepts_frame_values(self, value: object, expected: bool) -> None:
        assert from_records([_entry_record(decoy=value)])[0].decoy is expected

    def test_every_record_key_is_handled(self) -> None:
        # A column added to RECORD_KEYS without a home in from_records would be ignored,
        # which is the silent data loss this function exists to prevent.
        from pefftacular._records import _IDENTITY_KEYS, _SINGLE_KEYS, _VERBATIM_KEYS

        handled = {column for column, _ in _SINGLE_KEYS} | set(_IDENTITY_KEYS) | set(_VERBATIM_KEYS)
        assert handled == set(RECORD_KEYS)

    def test_the_peff_key_names_match_the_library(self) -> None:
        # _PEFF_KEYS above is this module's own copy; it must not drift from the library's.
        from pefftacular._records import _SINGLE_KEYS

        assert dict(_SINGLE_KEYS) == _PEFF_KEYS


def test_from_records_after_a_pandas_round_trip() -> None:
    pd = pytest.importorskip("pandas")
    header, entries = _read(FIXTURES / "custom_keys.peff")
    frame = pd.DataFrame(to_records(entries))
    assert from_records(frame.to_dict("records"), header) == entries
