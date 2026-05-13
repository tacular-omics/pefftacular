"""Frozen dataclass models for PEFF file structures."""

from dataclasses import dataclass, field
from datetime import date, time

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

    prefix: str | None = None
    db_name: str | None = None
    db_description: str | None = None
    db_version: str | None = None
    db_date: str | None = None
    db_sources: tuple[str, ...] = ()
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
    """Disulfide bond between two cysteine residues."""

    positions: tuple[int | str, ...]
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

    key_name: str
    fields: dict[str, CustomFieldValue] = field(default_factory=dict)
    raw: str = ""


# ---------------------------------------------------------------------------
# Sequence entry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SequenceEntry:
    """A single sequence entry in a PEFF file."""

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
