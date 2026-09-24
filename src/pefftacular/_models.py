"""Frozen dataclass models for PEFF file structures."""

from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import date, time
from typing import Literal, overload

CustomFieldValue = str | int | float | bool | date | time

# ---------------------------------------------------------------------------
# Header types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CustomKeyDef:
    """Definition of a custom key in a PEFF database header."""

    key_name: str
    description: str
    concept_curie: str | None = None
    regexp: str | None = None
    field_names: tuple[str, ...] = ()
    field_types: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class OptionalTagDef:
    """Definition of an optional tag used in entry annotations."""

    tag: str
    description: str


@dataclass(frozen=True, slots=True)
class DatabaseHeader:
    """Metadata block for a single database within a PEFF file."""

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    prefix: str | None = None
    db_name: str | None = None
    db_description: str | None = None
    db_version: str | None = None
    db_date: str | None = None
    db_sources: tuple[str, ...] = ()
    general_comments: tuple[str, ...] = ()
    number_of_entries: int | None = None
    sequence_type: str | None = None
    decoy: bool | None = None
    conversion: str | None = None
    has_annotation_identifiers: bool = False
    proteoform_db: bool = False
    custom_key_defs: tuple[CustomKeyDef, ...] = ()
    optional_tag_defs: tuple[OptionalTagDef, ...] = ()
    extra: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class FileHeader:
    """Top-level header for a PEFF file."""

    # Explicitly unhashable: its DatabaseHeaders hold an ``extra`` dict.
    __hash__ = None  # type: ignore[assignment]

    peff_version: str
    general_comments: tuple[str, ...] = ()
    databases: tuple[DatabaseHeader, ...] = ()


# ---------------------------------------------------------------------------
# Entry annotation types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VariantSimple:
    """Single amino-acid substitution."""

    position: int | str
    new_amino_acid: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class VariantComplex:
    """Multi-residue variant (insertion, deletion, or substitution)."""

    start_pos: int | str
    end_pos: int | str
    new_sequence: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModResUnimod:
    """Modification annotated with a UNIMOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModResPsi:
    """Modification annotated with a PSI-MOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class ModRes:
    """Modification with a generic accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class Processed:
    """Processed molecule annotation (signal peptide, mature protein, etc.)."""

    start_pos: int | str
    end_pos: int | str
    accession: str
    name: str
    tag: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class DisulfideBond:
    """Disulfide bond linking two previously-declared modification annotations.

    Per spec section 3.4.2, a ``\\DisulfideBond`` references two prior
    ``ModResPsi`` (half-cystine) entries *by their annotation IDs* — not by
    sequence position. ``annot_id_refs`` therefore holds those referenced IDs.
    """

    annot_id_refs: tuple[int | str, ...]
    description: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class SequenceRange:
    """A contiguous start–end range within a sequence (1-based, inclusive)."""

    start: int | str
    end: int | str


@dataclass(frozen=True, slots=True)
class Proteoform:
    """A specific proteoform — a combination of processing events and modifications."""

    proteoform_id: str
    ranges: tuple[SequenceRange, ...]
    annot_id_refs: tuple[int, ...]
    name: str | None = None
    annot_id: int | None = None


@dataclass(frozen=True, slots=True)
class CustomKeyValue:
    """A parsed value for a header-declared custom key on an entry."""

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    key_name: str
    fields: dict[str, CustomFieldValue] = field(default_factory=dict)
    raw: str = ""


# ---------------------------------------------------------------------------
# Sequence entry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SequenceEntry:
    """A single sequence entry in a PEFF file.

    Frozen and compared by value (``==``), but **not hashable**: ``custom_values`` and
    ``extra`` are dicts, so ``hash(entry)`` raises :class:`TypeError`. Key sets or dicts
    by ``(entry.prefix, entry.db_unique_id)`` instead.

    ``id`` holds the ``\\ID=`` key (spec §3.3.4); the name shadows the builtin but is kept
    for API stability.
    """

    # Explicitly unhashable (a generated hash would fail on the dict fields anyway).
    __hash__ = None  # type: ignore[assignment]

    prefix: str
    db_unique_id: str
    sequence: str
    id: str | None = None
    db_unique_id_key: str | None = None
    pname: str | None = None
    gname: str | None = None
    ncbi_tax_id: int | None = None
    tax_name: str | None = None
    length: int | None = None
    sv: int | None = None
    ev: int | None = None
    pe: int | None = None
    decoy: bool | None = None
    comment: str | None = None
    variant_simple: tuple[VariantSimple, ...] = ()
    variant_complex: tuple[VariantComplex, ...] = ()
    mod_res_unimod: tuple[ModResUnimod, ...] = ()
    mod_res_psi: tuple[ModResPsi, ...] = ()
    mod_res: tuple[ModRes, ...] = ()
    processed: tuple[Processed, ...] = ()
    disulfide_bond: tuple[DisulfideBond, ...] = ()
    proteoform: tuple[Proteoform, ...] = ()
    custom_values: dict[str, tuple[CustomKeyValue, ...]] = field(default_factory=dict)
    extra: dict[str, str] = field(default_factory=dict)

    # -- conversions (logic in _convert.py; models stay plain data) ---------------------

    @classmethod
    def from_fasta(cls, header: str, sequence: str, *, prefix: str | None = None) -> "SequenceEntry":
        """Build an entry from a plain FASTA header and sequence.

        pefftacular has no dependencies, so this takes strings rather than a
        ``fastatacular.SequenceEntry``; pass ``(e.raw_header, e.sequence)`` for one of those.
        UniProt-style headers are split into fields: ``sp|P12345|NAME_HUMAN Name OS=Homo
        sapiens OX=9606 GN=ABC PE=1 SV=2`` gives ``prefix="sp"``, ``db_unique_id="P12345"``,
        ``id="NAME_HUMAN"``, ``pname``, ``tax_name``, ``ncbi_tax_id``, ``gname``, ``pe``,
        ``sv``; other ``KEY=value`` pairs go to ``extra``. ``length`` is set from the sequence.

        Args:
            header: The header line, with or without the leading ``>``.
            sequence: Residues; whitespace is removed.
            prefix: PEFF database prefix. Default: the ``db`` of a ``db:ACC`` (PEFF style,
                checked first) or ``db|ACC...`` identifier. Required for other identifiers.

        Raises:
            PeffError: The header is empty, or no prefix is given or derivable.
        """
        from pefftacular._convert import entry_from_fasta

        return entry_from_fasta(cls, header, sequence, prefix=prefix)

    def to_fasta(self) -> tuple[str, str]:
        """Return ``(header, sequence)`` for a plain FASTA record (header without ``>``).

        The header is ``prefix|db_unique_id[|id] [pname] [OS=] [OX=] [GN=] [PE=] [SV=]``
        followed by ``extra`` as ``KEY=value``. When ``db_unique_id`` contains ``|`` the
        identifier is the PEFF form ``prefix:db_unique_id`` instead, and ``id`` is dropped.
        Annotations (variants, modifications, processing, proteoforms), ``comment``, ``ev``
        and ``decoy`` have no FASTA form and are dropped.
        ``SequenceEntry.from_fasta(*entry.to_fasta())`` restores the fields above as long as
        ``pname`` and the values contain no `` KEY=`` text.
        """
        from pefftacular._convert import entry_to_fasta

        return entry_to_fasta(self)

    @overload
    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = ...,
        variants: Iterable[VariantSimple] = ...,
        errors: Literal["raise"] = ...,
    ) -> str: ...

    @overload
    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = ...,
        variants: Iterable[VariantSimple] = ...,
        errors: Literal["skip"],
    ) -> str | None: ...

    def to_proforma(
        self,
        *,
        mods: Literal["psimod", "unimod"] = "psimod",
        variants: Iterable[VariantSimple] = (),
        errors: Literal["raise", "skip"] = "raise",
    ) -> str | None:
        """Render the sequence with its modifications as a ProForma 2.0 string.

        ``mods="psimod"`` writes ``\\ModResPsi`` sites (``S[MOD:00046]``) and
        ``"unimod"`` writes ``\\ModResUnimod`` sites (``S[UNIMOD:21]``). ``\\ModRes``
        sites are included when their accession is from the same vocabulary; others are
        left out. An empty accession is written by name (``[M:name]`` / ``[U:name]``). Every
        listed site is modified at once. Unknown positions (``?``) become a ProForma
        unknown-position prefix (``[MOD:00046]^2?SEQ``); they are dropped when a ``*``
        variant truncates the sequence, since they may lie in the removed part. A site
        listed in both ``\\ModResPsi``/``\\ModResUnimod`` and ``\\ModRes`` is written once.
        PEFF cannot tell a terminal modification from one on the terminal residue, so all
        are written on the residue.

        Real files contain entries whose sites lie past the end of the sequence (12 of
        the 20,431 entries of the neXtProt human PEFF). Converting a whole file, pass
        ``errors="skip"`` to get ``None`` for those entries instead of an exception::

            forms = [p for e in entries if (p := e.to_proforma(errors="skip")) is not None]

        Args:
            mods: Which vocabulary to render.
            variants: ``VariantSimple`` substitutions to apply, usually a subset of
                ``self.variant_simple``. A modification on a substituted residue is
                dropped (spec section 3.3.10: a modified variant needs its own entry).
                ``*`` (stop) truncates the sequence before that position.
            errors: ``"raise"`` (default) raises :class:`PeffError` for an entry that
                cannot be written; ``"skip"`` returns ``None`` for it instead.

        Returns:
            The ProForma string, or ``None`` with ``errors="skip"`` when the entry
            cannot be written.

        Raises:
            PeffError: Unknown ``mods`` or ``errors`` (always). With ``errors="raise"``
                also a position outside the sequence, a non-numeric position other than
                ``?``, two different substitutions at one position, or a square bracket
                in a modification; the message starts with ``prefix:db_unique_id``.
        """
        from pefftacular._convert import entry_to_proforma

        return entry_to_proforma(self, mods=mods, variants=variants, errors=errors)
