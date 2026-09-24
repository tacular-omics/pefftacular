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
        if acc.isdigit():
            acc = f"{cv}:{acc}"
        tag = acc
    else:
        tag = f"{name_prefix}:{mod.name}"
    if "[" in tag or "]" in tag:
        raise PeffError(
            f"{name}: modification {tag!r} contains a square bracket and cannot be written as ProForma",
        )
    return tag


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
        if len(new) != 1 or not new.isalpha():
            raise PeffError(f"{name}: VariantSimple at {pos}: new amino acid {new!r} is not a single letter or '*'")
        if substituted.get(pos, new) != new:
            raise PeffError(f"{name}: two different VariantSimple substitutions at position {pos}")
        substituted[pos] = new
    for pos, new in substituted.items():
        seq[pos - 1] = new

    own = entry.mod_res_psi if mods == "psimod" else entry.mod_res_unimod
    generic = tuple(m for m in entry.mod_res if m.accession.strip().upper().startswith(cv + ":"))
    at: dict[int, list[str]] = {}
    unknown: dict[str, int] = {}
    for source in (own, generic):
        # The same site may be listed in both \ModResPsi/\ModResUnimod and \ModRes:
        # write it once. Unknown sites count per list; the larger count wins.
        source_unknown: dict[str, int] = {}
        for mod in source:
            tag = _mod_tag(mod, cv, name_prefix, name)
            for p in mod.positions:
                if p == "?":
                    source_unknown[tag] = source_unknown.get(tag, 0) + 1
                    continue
                pos = _check_position(p, n, type(mod).__name__, name)
                if pos in substituted or pos > end:
                    continue  # the modified residue is replaced or truncated by a variant
                tags = at.setdefault(pos, [])
                if tag not in tags:
                    tags.append(tag)
        for tag, count in source_unknown.items():
            unknown[tag] = max(unknown.get(tag, 0), count)
    if end < n:
        unknown = {}  # an unknown site may lie in the truncated part: drop it

    head = "".join(f"[{tag}]" + (f"^{count}" if count > 1 else "") for tag, count in unknown.items())
    if head:
        head += "?"
    body = "".join(aa + "".join(f"[{tag}]" for tag in at.get(i, ())) for i, aa in enumerate(seq[:end], 1))
    return head + body


def _check_position(pos: int | str, n: int, what: str, name: str) -> int:
    if not isinstance(pos, int) or not 1 <= pos <= n:
        raise PeffError(f"{name}: {what} position {pos!r} is not a residue position (1..{n})")
    return pos
