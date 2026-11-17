r"""Tests for applying ``\VariantSimple`` and ``\VariantComplex`` to a sequence.

The structure follows the validation plan of the design: the spec's own §3.3.9 table,
position-map invariants over every shipped fixture, annotation remapping, agreement with
``to_proforma``, a write round trip, and the form counts.
"""

from __future__ import annotations

import dataclasses
import glob
import io
import math
import warnings
from pathlib import Path

import pytest

from pefftacular import (
    PeffError,
    PeffParseError,
    PeffWarning,
    SequenceEntry,
    VariantComplex,
    VariantForm,
    VariantSimple,
    read_peff,
    write_peff,
)

FIXTURES = Path(__file__).parent / "fixtures"
INSULIN = FIXTURES / "PEFF_AnnotID_Insulin_Valid.peff"

_HEADER = (
    "# PEFF 1.0\n# //\n# DbName=t\n# Prefix=t\n# DbVersion=1\n# DbSource=x\n# NumberOfEntries=1\n# SequenceType=AA\n"
)
SEQ = "ACDEFGHIKL"  # 10 residues


def _entry(description: str, *, sequence: str = SEQ, header_extra: str = "") -> SequenceEntry:
    data = f"{_HEADER}{header_extra}# //\n>t:X1 {description}\n{sequence}\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PeffWarning)
        return read_peff(io.StringIO(data))[1][0]


def _quiet(entry: SequenceEntry, variants) -> VariantForm:
    """Apply variants and fail the test if anything warns."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", PeffWarning)
        return entry.apply_variants(variants)


def _fixture_entries():
    """Every entry of every readable shipped fixture."""
    paths = sorted(glob.glob(str(FIXTURES / "*.peff")) + glob.glob(str(FIXTURES / "spec" / "*.peff")))
    for path in paths:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PeffWarning)
            try:
                _, entries = read_peff(path)
            except PeffParseError:
                continue  # the deliberately invalid fixtures
        for entry in entries:
            yield Path(path).name, entry


# ---------------------------------------------------------------------------
# 1. The spec's own examples (section 3.3.9 table, section 3.3.8 for '*')
# ---------------------------------------------------------------------------


class TestSpecExamples:
    """The §3.3.9 table, rewritten onto a 10-mer (positions 4..6 for 100..102)."""

    @pytest.mark.parametrize(
        ("variant", "expected"),
        [
            # "(100|100|)": a single amino acid deletion, with no character standing for it.
            (VariantComplex(4, 4, ""), "ACDFGHIKL"),
            # "(100|100||10kexomes)": the same, tagged.
            (VariantComplex(4, 4, "", tag="10kexomes"), "ACDFGHIKL"),
            # "(100|102|)": a 3-AA deletion starting at position 100.
            (VariantComplex(4, 6, ""), "ACDHIKL"),
            # "(100|100|APT)": replacement of the residue at 100 by APT; for an insertion the
            # inserted residues come after the residue defined in the sequence line.
            (VariantComplex(6, 6, "APT"), "ACDEFAPTHIKL"),
            # "(100|102|KPA)": a 3-AA substitution as a cassette.
            (VariantComplex(4, 6, "KPA"), "ACDKPAHIKL"),
            # "(100|101|P)": both residues removed and replaced with a single P.
            (VariantComplex(4, 5, "P"), "ACDPGHIKL"),
        ],
    )
    def test_variant_complex_table(self, variant: VariantComplex, expected: str) -> None:
        entry = _entry(r"\PName=x")
        assert _quiet(entry, [variant]).sequence == expected

    def test_illegal_single_residue_complex_is_applied_with_a_warning(self) -> None:
        # "(100|100|A) ILLEGAL": this must be a VariantSimple (§3.3.9). Reading is
        # permissive, so it is applied and reported rather than refused.
        entry = _entry(r"\PName=x")
        with pytest.warns(PeffWarning, match="should be a VariantSimple"):
            form = entry.apply_variants([VariantComplex(4, 4, "W")])
        assert form.sequence == "ACDWFGHIKL"

    def test_variant_simple_substitutes_one_residue(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)")
        form = _quiet(entry, entry.variant_simple)
        assert form.sequence == "ACDWFGHIKL"
        assert len(form.sequence) == len(SEQ)  # §3.3.8: MUST NOT extend the protein

    def test_stop_truncates_before_the_position(self) -> None:
        # §3.3.8: '*' is a nonsense mutation over which a peptide MUST NOT span.
        entry = _entry(r"\VariantSimple=(4|*)")
        form = _quiet(entry, entry.variant_simple)
        assert (form.sequence, form.truncated_at) == ("ACD", 4)
        assert form.position_map == (None, 1, 2, 3, None, None, None, None, None, None, None)

    def test_non_single_letter_substitution_raises(self) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match="is not a single letter A-Z or '\\*'"):
            entry.apply_variants([VariantSimple(4, "AC")])


# ---------------------------------------------------------------------------
# Positions, bounds and overlaps
# ---------------------------------------------------------------------------


class TestPositionRules:
    def test_position_outside_the_sequence_raises(self) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match=r"VariantSimple position 11 is outside 1\.\.10"):
            entry.apply_variants([VariantSimple(11, "W")])

    def test_complex_range_outside_the_sequence_raises(self) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match=r"VariantComplex range 9-11 is outside 1\.\.10"):
            entry.apply_variants([VariantComplex(9, 11, "W")])

    def test_reversed_complex_range_raises(self) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match="VariantComplex range 6-4"):
            entry.apply_variants([VariantComplex(6, 4, "W")])

    def test_unknown_position_cannot_locate_a_variant(self) -> None:
        # The reader accepts '?' in any position field; a variant cannot be applied there.
        entry = _entry(r"\VariantSimple=(?|W)")
        assert entry.variant_simple[0].position == "?"
        with pytest.raises(PeffError, match="position '\\?' is not a number"):
            entry.apply_variants(entry.variant_simple)

    def test_overlapping_variants_raise(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)(4|Y)")
        with pytest.raises(PeffError, match="both change position 4"):
            entry.apply_variants(entry.variant_simple)

    def test_overlap_between_simple_and_complex_raises(self) -> None:
        entry = _entry(r"\VariantSimple=(5|W) \VariantComplex=(4|6|KK)")
        with pytest.raises(PeffError, match="both change position 5"):
            entry.apply_variants([*entry.variant_simple, *entry.variant_complex])

    def test_adjacent_variants_compose(self) -> None:
        entry = _entry(r"\VariantSimple=(2|W) \VariantComplex=(5|6|)")
        form = _quiet(entry, [*entry.variant_simple, *entry.variant_complex])
        assert form.sequence == "AWDEHIKL"
        assert len(form.applied) == 2

    def test_errors_skip_returns_none(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)(4|Y)")
        assert entry.apply_variants(entry.variant_simple, errors="skip") is None

    def test_unknown_errors_mode_raises(self) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match="unknown errors mode"):
            entry.apply_variants([], errors="nope")  # type: ignore[arg-type]

    def test_no_variants_is_the_reference_form(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)")
        form = entry.apply_variants()
        assert form.is_reference
        assert form.entry is entry
        assert form.sequence == SEQ


# ---------------------------------------------------------------------------
# 2. Position-map invariants, over every shipped fixture
# ---------------------------------------------------------------------------


class TestPositionMap:
    def test_reference_map_is_the_identity(self) -> None:
        entry = _entry(r"\PName=x")
        form = entry.apply_variants()
        assert form.position_map == (None, *range(1, len(SEQ) + 1))

    def test_map_is_consistent_for_every_single_variant_form(self) -> None:
        checked = 0
        for name, entry in _fixture_entries():
            for variant in (*entry.variant_simple, *entry.variant_complex):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PeffWarning)
                    form = entry.apply_variants([variant], errors="skip")
                if form is None:
                    continue
                checked += 1
                where = f"{name} {entry.db_unique_id} {variant}"
                mapping = form.position_map
                assert len(mapping) == len(entry.sequence) + 1, where
                assert mapping[0] is None, where

                defined = [value for value in mapping[1:] if value is not None]
                assert defined == sorted(defined), where  # order is preserved
                assert len(set(defined)) == len(defined), where  # and nothing doubles up
                assert all(1 <= value <= len(form.sequence) for value in defined), where

                if form.truncated_at is not None:
                    assert len(form.sequence) == form.truncated_at - 1, where
                    assert all(value is None for value in mapping[form.truncated_at :]), where
                    continue

                if isinstance(variant, VariantSimple):
                    assert len(form.sequence) == len(entry.sequence), where
                    assert None not in mapping[1:], where
                else:
                    start, end = int(variant.start_pos), int(variant.end_pos)
                    delta = len(variant.new_sequence) - (end - start + 1)
                    assert len(form.sequence) == len(entry.sequence) + delta, where
                    # None exactly for the replaced span.
                    assert [pos for pos, value in enumerate(mapping) if value is None and pos] == list(
                        range(start, end + 1)
                    ), where
        assert checked > 800, f"only {checked} forms checked"


# ---------------------------------------------------------------------------
# 3. Annotation remapping
# ---------------------------------------------------------------------------


class TestAnnotationRemapping:
    MODS = (
        r"\ModResPsi=(2|MOD:00046|O-phospho-L-serine)(5|MOD:00046|O-phospho-L-serine)(9|MOD:00046|O-phospho-L-serine)"
    )

    def test_sites_before_inside_and_after_a_deletion(self) -> None:
        entry = _entry(f"{self.MODS} \\VariantComplex=(4|6|)")
        with pytest.warns(PeffWarning, match="the residue was removed"):
            form = entry.apply_variants(entry.variant_complex)
        positions = [p for mod in form.entry.mod_res_psi for p in mod.positions]
        assert positions == [2, 6]  # 2 unchanged, 5 removed, 9 shifted by -3

    def test_sites_after_an_insertion_shift_right(self) -> None:
        # Nothing is annotated at position 6, so nothing is dropped and nothing warns.
        entry = _entry(f"{self.MODS} \\VariantComplex=(6|6|APT)")
        form = _quiet(entry, entry.variant_complex)
        positions = [p for mod in form.entry.mod_res_psi for p in mod.positions]
        assert positions == [2, 5, 11]  # 9 + 2 inserted residues

    def test_site_on_a_substituted_residue_is_dropped(self) -> None:
        # §3.3.10: a modified variant needs its own entry, so the site cannot carry over.
        entry = _entry(f"{self.MODS} \\VariantSimple=(5|W)")
        with pytest.warns(PeffWarning, match="the residue was substituted"):
            form = entry.apply_variants(entry.variant_simple)
        positions = [p for mod in form.entry.mod_res_psi for p in mod.positions]
        assert positions == [2, 9]

    def test_unknown_position_is_kept_as_written(self) -> None:
        entry = _entry(r"\ModResPsi=(?|MOD:00046|O-phospho-L-serine) \VariantComplex=(4|6|)")
        form = _quiet(entry, entry.variant_complex)
        assert form.entry.mod_res_psi[0].positions == ("?",)

    def test_modification_with_no_remaining_site_is_dropped(self) -> None:
        entry = _entry(r"\ModResPsi=(5|MOD:00046|O-phospho-L-serine) \VariantComplex=(4|6|)")
        with pytest.warns(PeffWarning, match="the residue was removed"):
            form = entry.apply_variants(entry.variant_complex)
        assert form.entry.mod_res_psi == ()

    def test_processed_range_is_shifted(self) -> None:
        entry = _entry(r"\Processed=(7|10|PEFF:0001021|mature protein) \VariantComplex=(2|3|)")
        form = _quiet(entry, entry.variant_complex)
        assert (form.entry.processed[0].start_pos, form.entry.processed[0].end_pos) == (5, 8)

    def test_processed_range_losing_a_bound_is_dropped_not_clipped(self) -> None:
        entry = _entry(r"\Processed=(5|10|PEFF:0001021|mature protein) \VariantComplex=(4|6|)")
        with pytest.warns(PeffWarning, match="Processed .* a bound was removed"):
            form = entry.apply_variants(entry.variant_complex)
        assert form.entry.processed == ()

    def test_pending_variants_are_remapped_and_the_applied_one_removed(self) -> None:
        entry = _entry(r"\VariantSimple=(2|W)(9|Y) \VariantComplex=(4|6|)")
        form = _quiet(entry, entry.variant_complex)
        assert [v.position for v in form.entry.variant_simple] == [2, 6]
        assert form.entry.variant_complex == ()
        assert form.applied == entry.variant_complex

    def test_length_is_recomputed_only_when_the_entry_has_one(self) -> None:
        with_length = _entry(r"\Length=10 \VariantComplex=(4|6|)")
        assert _quiet(with_length, with_length.variant_complex).entry.length == 7
        without = _entry(r"\VariantComplex=(4|6|)")
        assert _quiet(without, without.variant_complex).entry.length is None

    def test_disulfide_bond_refs_are_ids_not_positions(self) -> None:
        entry = _entry(
            r"\ModResPsi=(1:2|MOD:00798|half cystine)(2:9|MOD:00798|half cystine) "
            r"\DisulfideBond=(3:1,2|between chains) \VariantComplex=(4|6|)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        form = _quiet(entry, entry.variant_complex)
        assert form.entry.disulfide_bond[0].annot_id_refs == (1, 2)
        assert [p for mod in form.entry.mod_res_psi for p in mod.positions] == [2, 6]

    def test_annotation_whose_reference_is_gone_is_dropped(self) -> None:
        entry = _entry(
            r"\ModResPsi=(1:5|MOD:00798|half cystine)(2:9|MOD:00798|half cystine) "
            r"\DisulfideBond=(3:1,2|between chains) \VariantComplex=(4|6|)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        with pytest.warns(PeffWarning, match="the residue was removed"):
            form = entry.apply_variants(entry.variant_complex)
        assert form.entry.mod_res_psi[0].positions == (6,)  # id 1 is gone, id 2 shifted
        assert form.entry.disulfide_bond == ()

    def test_proteoform_referencing_an_applied_variant_is_dropped(self) -> None:
        entry = _entry(
            r"\VariantSimple=(1:4|W) \Proteoform=(pf1|1-10|1|with the variant)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        form = _quiet(entry, entry.variant_simple)
        assert form.sequence == "ACDWFGHIKL"
        assert form.entry.proteoform == ()

    def test_proteoform_not_referencing_the_variant_is_kept_and_remapped(self) -> None:
        entry = _entry(
            r"\ModResPsi=(1:9|MOD:00046|O-phospho-L-serine) \VariantComplex=(4|6|) "
            r"\Proteoform=(pf1|7-10|1|tail)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        form = _quiet(entry, entry.variant_complex)
        (proteoform,) = form.entry.proteoform
        assert (proteoform.ranges[0].start, proteoform.ranges[0].end) == (4, 7)
        assert proteoform.annot_id_refs == (1,)

    def test_site_already_past_the_sequence_end_is_dropped(self) -> None:
        # The reader warns about such a site but keeps it; it cannot be remapped.
        entry = _entry(r"\ModResPsi=(99|MOD:00046|O-phospho-L-serine) \VariantComplex=(4|6|)")
        with pytest.warns(PeffWarning, match="the residue was removed"):
            form = entry.apply_variants(entry.variant_complex)
        assert form.entry.mod_res_psi == ()

    def test_dropped_processed_takes_its_references_with_it(self) -> None:
        entry = _entry(
            r"\Processed=(1:5|10|PEFF:0001021|mature protein) \VariantComplex=(4|6|) "
            r"\Proteoform=(pf1|1-7|1|the mature form)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        with pytest.warns(PeffWarning, match="Processed .* a bound was removed"):
            form = entry.apply_variants(entry.variant_complex)
        assert form.entry.processed == ()
        assert form.entry.proteoform == ()  # its only reference is gone

    def test_other_fields_and_custom_values_survive(self) -> None:
        entry = _entry(r"\PName=Insulin \GName=INS \NcbiTaxId=9606 \Unknown=x \VariantComplex=(4|6|)")
        form = _quiet(entry, entry.variant_complex)
        assert (form.entry.pname, form.entry.gname, form.entry.ncbi_tax_id) == ("Insulin", "INS", 9606)
        assert form.entry.extra == {"Unknown": "x"}
        assert (form.entry.prefix, form.entry.db_unique_id) == ("t", "X1")


# ---------------------------------------------------------------------------
# 4. Agreement with to_proforma, which already applies VariantSimple
# ---------------------------------------------------------------------------


def test_apply_variants_agrees_with_to_proforma_on_every_fixture() -> None:
    """``to_proforma(variants=[v])`` and ``apply_variants([v])`` must render the same."""
    compared = 0
    for name, entry in _fixture_entries():
        for variant in entry.variant_simple:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PeffWarning)
                left = entry.to_proforma(variants=[variant], errors="skip")
                form = entry.apply_variants([variant], errors="skip")
            right = form.entry.to_proforma(errors="skip") if form is not None else None
            assert left == right, f"{name} {entry.db_unique_id} {variant}"
            compared += 1
    assert compared > 800, f"only {compared} pairs compared"


# ---------------------------------------------------------------------------
# 5. A form writes and reads back
# ---------------------------------------------------------------------------


def test_form_round_trips_through_write_peff() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", PeffWarning)
        header, entries = read_peff(FIXTURES / "PEFF_Tiny_Valid.peff")
    entry = next(e for e in entries if e.variant_simple)
    form = entry.apply_variants(entry.variant_simple[:1])

    # Write just this entry, so the header must declare just its database.
    database = next(db for db in header.databases if db.prefix == entry.prefix)
    single = dataclasses.replace(header, databases=(dataclasses.replace(database, number_of_entries=1),))
    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("error", PeffWarning)
        write_peff(single, [form.entry], buf)
        _, again = read_peff(io.StringIO(buf.getvalue()))
    assert again == [form.entry]


# ---------------------------------------------------------------------------
# 6. How many forms, and the limit
# ---------------------------------------------------------------------------


class TestFormCounts:
    """The numbers measured on the spec's insulin example (§3.4.2)."""

    @staticmethod
    def _insulin() -> SequenceEntry:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PeffWarning)
            return read_peff(INSULIN)[1][0]

    def test_the_fixture_still_has_the_annotations_the_counts_assume(self) -> None:
        entry = self._insulin()
        assert len(entry.variant_simple) == 70
        assert entry.variant_complex == ()
        assert len(entry.proteoform) == 11

    def test_default_is_the_reference_plus_one_form_per_variant(self) -> None:
        entry = self._insulin()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PeffWarning)
            forms = list(entry.variant_forms())
        assert len(forms) == 71
        assert forms[0].is_reference
        assert all(len(form.applied) == 1 for form in forms[1:])

    def test_pairs_count_2486_and_yield_2466_because_20_pairs_overlap(self) -> None:
        entry = self._insulin()
        assert 1 + 70 + math.comb(70, 2) == 2486  # what the limit is checked against
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", PeffWarning)
            forms = list(entry.variant_forms(max_combination=2, limit=2486))
        skipped = [w for w in caught if "combination skipped" in str(w.message)]
        # The file lists alternative substitutions at one position, which cannot combine.
        assert len(skipped) == 20
        assert len(forms) == 2486 - 20 == 2466

    def test_proteoforms_yields_the_combinations_the_file_declares(self) -> None:
        entry = self._insulin()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PeffWarning)
            forms = list(entry.variant_forms(proteoforms=True))
        assert len(forms) == 11

    def test_a_proteoform_of_contradictory_variants_is_skipped(self) -> None:
        entry = _entry(
            r"\VariantSimple=(1:4|W)(2:4|Y) \Proteoform=(pf1|1-10|1,2|both at once)",
            header_extra="# HasAnnotationIdentifiers=true\n",
        )
        with pytest.warns(PeffWarning, match="proteoform pf1 skipped: .* both change position 4"):
            forms = list(entry.variant_forms(proteoforms=True))
        assert forms == []

    def test_limit_raises_before_any_sequence_is_built(self) -> None:
        entry = self._insulin()
        iterator = entry.variant_forms(max_combination=2)  # default limit 10,000 < 2486? no: 2486 < 10000
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PeffWarning)
            assert len(list(iterator)) == 2466
        with pytest.raises(PeffError, match="over the limit of 100"):
            next(entry.variant_forms(max_combination=2, limit=100))
        with pytest.raises(PeffError, match="70 variants give 71 forms"):
            next(entry.variant_forms(limit=70))

    def test_power_set_is_refused_rather_than_attempted(self) -> None:
        entry = self._insulin()
        with pytest.raises(PeffError, match="max_combination=70"):
            next(entry.variant_forms(max_combination=70))

    def test_negative_max_combination_raises(self) -> None:
        entry = self._insulin()
        with pytest.raises(PeffError, match="must not be negative"):
            next(entry.variant_forms(max_combination=-1))

    def test_include_reference_false_drops_the_first_form(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)(9|Y)")
        forms = list(entry.variant_forms(include_reference=False))
        assert len(forms) == 2
        assert not any(form.is_reference for form in forms)

    def test_entry_without_variants_yields_only_the_reference(self) -> None:
        entry = _entry(r"\PName=x")
        (form,) = list(entry.variant_forms())
        assert form.is_reference and form.sequence == SEQ

    def test_mandatory_applies_every_variant_in_one_form(self) -> None:
        # §3.4.1: with ProteoformDb=true the annotations MUST all be applied.
        entry = _entry(r"\VariantSimple=(2|W) \VariantComplex=(5|6|)", header_extra="# ProteoformDb=true\n")
        (form,) = list(entry.variant_forms(mandatory=True))
        assert form.sequence == "AWDEHIKL"
        assert len(form.applied) == 2

    def test_mandatory_on_contradictory_variants_raises(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)(4|Y)")
        with pytest.raises(PeffError, match="both change position 4"):
            next(entry.variant_forms(mandatory=True))

    def test_forms_are_lazy(self) -> None:
        entry = self._insulin()
        iterator = entry.variant_forms(max_combination=3, limit=10**6)
        assert next(iterator).is_reference  # never builds the 57,226 forms


# ---------------------------------------------------------------------------
# 7. What counts as applied, and what may follow a stop
# ---------------------------------------------------------------------------


class TestWhatCountsAsApplied:
    """The models are frozen values, so an equal variant is the same variant."""

    def test_an_equal_variant_is_the_one_applied(self) -> None:
        entry = _entry(r"\VariantSimple=(4|W)")
        rebuilt = VariantSimple(4, "W")
        assert rebuilt == entry.variant_simple[0] and rebuilt is not entry.variant_simple[0]
        form = _quiet(entry, [rebuilt])
        assert form.sequence == "ACDWFGHIKL"
        # It is in the sequence now, so it must not still be listed as pending.
        assert form.entry.variant_simple == ()

    def test_an_equal_complex_variant_is_the_one_applied(self) -> None:
        entry = _entry(r"\VariantComplex=(4|6|KPA)")
        form = _quiet(entry, [VariantComplex(4, 6, "KPA")])
        assert (form.sequence, form.entry.variant_complex) == ("ACDKPAHIKL", ())

    def test_a_variant_listed_twice_keeps_the_second_as_pending(self) -> None:
        # One application consumes one occurrence, so the duplicate survives.
        entry = _entry(r"\VariantSimple=(4|W)(4|W)")
        assert entry.variant_simple == (VariantSimple(4, "W"), VariantSimple(4, "W"))
        form = _quiet(entry, [entry.variant_simple[0]])
        assert form.entry.variant_simple == (VariantSimple(4, "W"),)


class TestNothingFollowsAStop:
    def test_a_variant_after_a_stop_raises(self) -> None:
        # The protein ends at the stop (spec section 3.3.8), so position 8 does not
        # exist in this form and cannot be substituted in it.
        entry = _entry(r"\VariantSimple=(4|*)(8|Q)")
        with pytest.raises(PeffError, match=r"VariantSimple\(8\|Q\) is after the stop VariantSimple\(4\|\*\) at 4"):
            entry.apply_variants(entry.variant_simple)

    def test_a_complex_variant_after_a_stop_raises(self) -> None:
        entry = _entry(r"\VariantSimple=(4|*) \VariantComplex=(7|8|KK)")
        with pytest.raises(PeffError, match="is after the stop"):
            entry.apply_variants([*entry.variant_simple, *entry.variant_complex])

    def test_a_stop_last_applies_every_variant_it_lists(self) -> None:
        entry = _entry(r"\VariantSimple=(2|Q)(4|*)")
        form = _quiet(entry, entry.variant_simple)
        assert (form.sequence, form.truncated_at) == ("AQD", 4)
        assert form.applied == entry.variant_simple  # every one of them is in the sequence

    def test_variant_forms_skips_the_impossible_combination(self) -> None:
        entry = _entry(r"\VariantSimple=(4|*)(8|Q)")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            forms = list(entry.variant_forms(max_combination=2))
        # The reference and the two single-variant forms; the pair is refused, not yielded.
        assert [form.sequence for form in forms] == [SEQ, "ACD", "ACDEFGHQKL"]
        skipped = [str(w.message) for w in caught if "combination skipped" in str(w.message)]
        assert len(skipped) == 1 and "is after the stop" in skipped[0]


class TestComplexReplacementIsValidated:
    """Spec section 3.3.9: residue letters, or nothing at all."""

    @pytest.mark.parametrize("new_sequence", ["Q*", "*", " Q ", "q1", "-", "[AEQ]P", "K P"])
    def test_a_replacement_that_is_not_residue_letters_raises(self, new_sequence: str) -> None:
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match="is not amino acid letters A-Z"):
            entry.apply_variants([VariantComplex(4, 5, new_sequence)])

    def test_the_asterisk_is_refused_rather_than_truncating(self) -> None:
        # Section 3.3.9 never gives the asterisk a meaning; the nonsense codon is a
        # VariantSimple (section 3.3.8). So this is refused instead of spliced in.
        entry = _entry(r"\PName=x")
        with pytest.raises(PeffError, match=r"no '\*' \(the stop codon is a VariantSimple"):
            entry.apply_variants([VariantComplex(4, 5, "Q*")])

    def test_an_empty_replacement_is_still_a_deletion(self) -> None:
        entry = _entry(r"\PName=x")
        assert _quiet(entry, [VariantComplex(4, 5, "")]).sequence == "ACDGHIKL"  # E and F removed


def test_variant_form_is_unhashable_like_the_entry() -> None:
    entry = _entry(r"\VariantSimple=(4|W)")
    form = entry.apply_variants(entry.variant_simple)
    assert isinstance(form, VariantForm)
    with pytest.raises(TypeError):
        hash(form)
