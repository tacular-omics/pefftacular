"""Conversions between PEFF entries and plain FASTA / ProForma.

pefftacular has no runtime dependencies, so FASTA is exchanged as a plain header string
and a sequence string rather than as ``fastatacular.SequenceEntry`` objects. The header
format matches UniProt and ``fastatacular``: ``db|ACCESSION|ENTRY_NAME protein name
OS=... OX=... GN=... PE=... SV=...``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import TYPE_CHECKING, Literal

from pefftacular.errors import PeffError

if TYPE_CHECKING:
    from pefftacular._models import ModRes, ModResPsi, ModResUnimod, SequenceEntry, VariantSimple

# Same rules as fastatacular: KEY=value pairs, the value running to the next " KEY=".
_KV_PATTERN = re.compile(r"(?:^|(?<=\s))(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<val>.*?)(?=\s+[A-Za-z_][A-Za-z0-9_]*=|$)")
_UNIPROT_ID = re.compile(r"^(?P<prefix>[^|\s]+)\|(?P<accession>[^|]+)\|(?P<entry_name>[^|\s]+)$")
# NCBI ``gi|<number>|<db>|<accession>|[<chain or name>|...]``: the accession is the 4th field,
# not the gi number. ``gi`` may carry a decoy/contaminant tag (``DECOY_gi``, ``rev_gi``,
# ``REV-2-gi``; same rule as fastatacular), kept on the db prefix as for other header styles
# (``DECOY_sp|...`` gives prefix ``DECOY_sp``).
_NCBI_GI_ID = re.compile(
    r"^(?P<tag>(?i:(?:DECOY|REVERSE|REV|CONTAM|CON)(?:_|-[0-9]+-)))?gi"
    r"\|[0-9]+\|(?P<prefix>[A-Za-z]+)\|(?P<accession>[^|\s]+)"
    r"(?:\|(?P<entry_name>[^|\s]*)(?:\|.*)?)?$"
)
_UNIPROT_DBS = frozenset({"sp", "tr"})
_PIPE_ID = re.compile(r"^(?P<prefix>[^|\s]+)\|(?P<accession>[^|\s]+)(?:\|.*)?$")
_KEY = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_INT_KEYS = frozenset({"OX", "PE", "SV"})
_STR_KEYS = frozenset({"OS", "GN"})
_WS = re.compile(r"\s+")


def entry_from_fasta(
    cls: type[SequenceEntry], header: str, sequence: str, *, prefix: str | None = None
) -> SequenceEntry:
    header = header.strip()
    if header.startswith(">"):
        header = header[1:].lstrip()
    if not header:
        raise PeffError("Empty FASTA header")
    identifier, *rest = header.split(None, 1)
    description = rest[0].strip() if rest else ""

    db: str | None
    entry_name: str | None = None
    peff_db, colon, peff_acc = identifier.partition(":")
    if colon and peff_db and peff_acc and "|" not in peff_db:  # PEFF style, e.g. nxp:NX_P01308-1
        db, accession = peff_db, peff_acc
    elif m := _NCBI_GI_ID.match(identifier):
        db, accession = (m["tag"] or "") + m["prefix"], m["accession"]
        if m["prefix"] in _UNIPROT_DBS:
            entry_name = m["entry_name"] or None
        elif m["prefix"].lower() == "pdb" and m["entry_name"]:  # NCBI writes PDB chains as pdb|1MBA|A -> 1MBA_A
            accession = f"{accession}_{m['entry_name']}"
    elif m := _UNIPROT_ID.match(identifier):
        db, accession, entry_name = m["prefix"], m["accession"], m["entry_name"]
    elif m := _PIPE_ID.match(identifier):
        db, accession = m["prefix"], m["accession"]
    else:
        db, accession = None, identifier
    db = prefix or db
    if not db:
        raise _no_prefix(identifier)

    strs: dict[str, str] = {}
    ints: dict[str, int] = {}
    extra: dict[str, str] = {}
    first = _KV_PATTERN.search(description)
    pname = (description[: first.start()] if first else description).strip() or None
    if first:
        for kv in _KV_PATTERN.finditer(description, first.start()):
            key, val = kv["key"], kv["val"].strip()
            if key in _STR_KEYS and key not in strs:
                strs[key] = val
            elif key in _INT_KEYS and key not in ints and _is_int(val):
                ints[key] = int(val)
            else:
                extra[key] = val
    seq = _WS.sub("", sequence)
    return cls(
        prefix=db,
        db_unique_id=accession,
        sequence=seq,
        id=entry_name,
        pname=pname,
        tax_name=strs.get("OS"),
        gname=strs.get("GN"),
        ncbi_tax_id=ints.get("OX"),
        pe=ints.get("PE"),
        sv=ints.get("SV"),
        length=len(seq),
        extra=extra,
    )


def _no_prefix(identifier: str) -> PeffError:
    err = PeffError(f"Cannot tell the PEFF prefix of FASTA identifier {identifier!r}")
    err.add_note("hint: pass prefix=, e.g. SequenceEntry.from_fasta(header, sequence, prefix='gen')")
    return err


def _is_int(val: str) -> bool:
    try:
        int(val)
    except ValueError:
        return False
    return True


def entry_to_fasta(entry: SequenceEntry) -> tuple[str, str]:
    if "|" in entry.db_unique_id or ":" in entry.prefix:
        # Pipe form would be ambiguous; use the PEFF identifier (entry.id has no place here).
        parts = [f"{entry.prefix}:{entry.db_unique_id}"]
    else:
        parts = [f"{entry.prefix}|{entry.db_unique_id}" + (f"|{entry.id}" if entry.id else "")]
    if entry.pname:
        parts.append(entry.pname)
    for key, value in (
        ("OS", entry.tax_name),
        ("OX", entry.ncbi_tax_id),
        ("GN", entry.gname),
        ("PE", entry.pe),
        ("SV", entry.sv),
    ):
        if value is not None:
            parts.append(f"{key}={value}")
    parts.extend(f"{key}={value}" for key, value in entry.extra.items() if _KEY.fullmatch(key) and "\n" not in value)
    return " ".join(parts), entry.sequence


ModVocabulary = Literal["psimod", "unimod"]
_VOCAB = {"psimod": ("MOD", "M"), "unimod": ("UNIMOD", "U")}


def _mod_tag(mod: ModResPsi | ModResUnimod | ModRes, cv: str, name_prefix: str, name: str) -> str:
    acc = mod.accession.strip()
    if acc:
        if acc.isdigit():  # bare number in \ModResPsi / \ModResUnimod: that list's CV
            acc = f"{cv}:{acc}"
        tag = _canonical_accession(acc)
    else:
        tag = f"{name_prefix}:{mod.name}"
    # Inside a ProForma tag "[" / "]" end it, "|" starts another tag (an alternative or
    # INFO:) and "#" a group label (ProForma 2.0 has no escape): [M:a#b] is "a" in group "b".
    bad = next((c for c in "[]|#" if c in tag), None)
    if bad is not None:
        what = "a square bracket" if bad in "[]" else repr(bad)
        raise PeffError(f"{name}: modification {tag!r} contains {what} and cannot be written as ProForma")
    return tag


# Numeric CV accessions written in canonical form, so equal accessions in different
# spellings ("MOD:46", "MOD: 00046", "mod:00046") become one tag: PSI-MOD ids are zero-padded
# to 5 digits, Unimod ids are not padded.
_NUMERIC_CV_FORMAT = {"MOD": "{:05d}", "UNIMOD": "{:d}"}


def _canonical_accession(acc: str) -> str:
    """``"mod: 46"`` -> ``"MOD:00046"``; other CVs keep their value, prefix upper-cased."""
    prefix, sep, value = acc.partition(":")
    prefix, value = prefix.strip().upper(), value.strip()
    if not sep:
        return acc
    fmt = _NUMERIC_CV_FORMAT.get(prefix)
    if fmt is not None and value.isascii() and value.isdigit():
        value = fmt.format(int(value))
    return f"{prefix}:{value}"


def _tag_key(tag: str, name_prefix: str) -> str:
    """Normalised identity of a tag: accessions ignore case, names (``M:``/``U:``) do not."""
    return tag if tag.startswith(name_prefix + ":") else tag.upper()


ProformaErrors = Literal["raise", "skip"]


def entry_to_proforma(
    entry: SequenceEntry,
    *,
    mods: ModVocabulary = "psimod",
    variants: Iterable[VariantSimple] = (),
    errors: ProformaErrors = "raise",
) -> str | None:
    if mods not in _VOCAB:
        raise PeffError(f"mods must be 'psimod' or 'unimod', not {mods!r}")
    if errors not in ("raise", "skip"):
        raise PeffError(f"errors must be 'raise' or 'skip', not {errors!r}")
    try:
        return _to_proforma(entry, mods, variants)
    except PeffError:
        if errors == "skip":
            return None
        raise


def _to_proforma(entry: SequenceEntry, mods: ModVocabulary, variants: Iterable[VariantSimple]) -> str:
    cv, name_prefix = _VOCAB[mods]
    name = f"{entry.prefix}:{entry.db_unique_id}"
    seq = list(entry.sequence)
    n = len(seq)

    substituted: dict[int, str] = {}
    end = n
    for v in variants:
        pos = _check_position(v.position, n, "VariantSimple", name)
        new = v.new_amino_acid
        if new == "*":  # stop codon: the protein ends before this position
            end = min(end, pos - 1)
            continue
        if len(new) != 1 or not ("A" <= new <= "Z"):
            # Checked here, not as a residue below, so the error names the variant.
            raise PeffError(f"{name}: VariantSimple at {pos}: new amino acid {new!r} is not a single letter A-Z or '*'")
        if substituted.get(pos, new) != new:
            raise PeffError(f"{name}: two different VariantSimple substitutions at position {pos}")
        substituted[pos] = new
    for pos, new in substituted.items():
        seq[pos - 1] = new

    own = entry.mod_res_psi if mods == "psimod" else entry.mod_res_unimod
    generic = tuple(m for m in entry.mod_res if _canonical_accession(m.accession.strip()).startswith(cv + ":"))
    at: dict[int, list[str]] = {}  # position -> written tag texts, in first-seen order
    unknown: dict[str, int] = {}
    tag_text: dict[str, str] = {}  # normalised key -> first spelling seen
    for source in (own, generic):
        # The same site may be listed in both \ModResPsi/\ModResUnimod and \ModRes, in any
        # spelling ("MOD:00046", "mod:00046", "MOD:46", "MOD: 00046"): write it once, or a
        # reader adds its mass twice. Numeric MOD/UNIMOD accessions are written in canonical
        # form; other tags are compared by their normalised key and the first spelling wins.
        # Unknown sites count per list; the larger count wins.
        source_unknown: dict[str, int] = {}
        for mod in source:
            tag = _mod_tag(mod, cv, name_prefix, name)
            text = tag_text.setdefault(_tag_key(tag, name_prefix), tag)
            for p in mod.positions:
                if p == "?":
                    source_unknown[text] = source_unknown.get(text, 0) + 1
                    continue
                pos = _check_position(p, n, type(mod).__name__, name)
                if pos in substituted or pos > end:
                    continue  # the modified residue is replaced or truncated by a variant
                tags = at.setdefault(pos, [])
                if text not in tags:
                    tags.append(text)
        for text, count in source_unknown.items():
            unknown[text] = max(unknown.get(text, 0), count)
    if end < n:
        unknown = {}  # an unknown site may lie in the truncated part: drop it

    for i, aa in enumerate(seq[:end], 1):
        if not ("A" <= aa <= "Z"):
            # PEFF also allows "*" (interruption) and "-" (gap); ProForma 2.0 residues are A-Z only.
            raise PeffError(f"{name}: residue {aa!r} at position {i} cannot be written as ProForma (A-Z only)")

    head = "".join(f"[{tag}]" + (f"^{count}" if count > 1 else "") for tag, count in unknown.items())
    if head:
        head += "?"
    body = "".join(aa + "".join(f"[{tag}]" for tag in at.get(i, ())) for i, aa in enumerate(seq[:end], 1))
    return head + body


def _check_position(pos: int | str, n: int, what: str, name: str) -> int:
    if not isinstance(pos, int) or not 1 <= pos <= n:
        raise PeffError(f"{name}: {what} position {pos!r} is not a residue position (1..{n})")
    return pos
