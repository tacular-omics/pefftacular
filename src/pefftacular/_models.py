"""Frozen dataclass models for PEFF file structures."""

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Header types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CustomKeyDef:
    """Definition of a custom key in a PEFF database header."""

    key_name: str
    description: str
    regexp: str | None = None
    field_names: tuple[str, ...] = ()
    field_types: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class DatabaseHeader:
    """Metadata block for a single database within a PEFF file."""

    prefix: str | None = None
    db_name: str | None = None
    db_version: str | None = None
    db_source: str | None = None
    number_of_entries: int | None = None
    sequence_type: str | None = None
    custom_key_defs: tuple[CustomKeyDef, ...] = ()
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


@dataclass(frozen=True, slots=True)
class VariantComplex:
    """Multi-residue variant (insertion, deletion, or substitution)."""

    start_pos: int | str
    end_pos: int | str
    new_sequence: str
    tag: str | None = None


@dataclass(frozen=True, slots=True)
class ModResUnimod:
    """Modification annotated with a UNIMOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None


@dataclass(frozen=True, slots=True)
class ModResPsi:
    """Modification annotated with a PSI-MOD accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None


@dataclass(frozen=True, slots=True)
class ModRes:
    """Modification with a generic accession."""

    positions: tuple[int | str, ...]
    accession: str
    name: str
    tag: str | None = None


@dataclass(frozen=True, slots=True)
class Processed:
    """Processed molecule annotation (signal peptide, mature protein, etc.)."""

    start_pos: int | str
    end_pos: int | str
    accession: str
    name: str
    tag: str | None = None


# ---------------------------------------------------------------------------
# Sequence entry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SequenceEntry:
    """A single sequence entry in a PEFF file."""

    prefix: str
    db_unique_id: str
    sequence: str
    pname: str | None = None
    gname: str | None = None
    ncbi_tax_id: int | None = None
    tax_name: str | None = None
    length: int | None = None
    sv: int | None = None
    ev: int | None = None
    pe: int | None = None
    decoy: bool | None = None
    variant_simple: tuple[VariantSimple, ...] = ()
    variant_complex: tuple[VariantComplex, ...] = ()
    mod_res_unimod: tuple[ModResUnimod, ...] = ()
    mod_res_psi: tuple[ModResPsi, ...] = ()
    mod_res: tuple[ModRes, ...] = ()
    processed: tuple[Processed, ...] = ()
    extra: dict[str, str] = field(default_factory=dict)
