"""Turn UniProt JSON entries into labelled entries and precursor regions.

Regions (1-based inclusive UniProt coordinates):
  SIG   'Signal' feature
  PRO   'Propeptide' ending before the first mature peptide (the one
        adjacent to it if several)
  MAT   'Peptide' features (mature conopeptides); 'Chain' features only when
        an entry has no 'Peptide' and the chain does not overlap the signal
        or a propeptide
  POST  'Propeptide' starting after the last mature peptide
  INTER propeptides between mature peptides, or extra N-terminal ones
        (recorded, not used for profiles)

Labels come from the SIMILARITY comment ("Belongs to the conotoxin C
superfamily. Consomatin family.") and are normalised by built-in rules, then
by the curated overrides file, which always wins.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

from .config import LabelsConfig, RegionsConfig
from .core import cysteine_framework  # noqa: F401  (re-exported; defined in core)
from .io import read_tsv

log = logging.getLogger(__name__)

REGIONS = ("SIG", "PRO", "MAT", "POST", "INTER")
PROFILE_REGIONS = ("SIG", "PRO", "MAT", "POST")

# Evidence classes, strongest first.
_EVIDENCE = [
    ("experimental", {"ECO:0000269", "ECO:0007744", "ECO:0000314", "ECO:0000315"}),
    ("curated", {"ECO:0000305", "ECO:0000303", "ECO:0000312"}),
    ("similarity", {"ECO:0000250", "ECO:0000266"}),
    ("predicted", {"ECO:0000255", "ECO:0000256", "ECO:0000259"}),
]
EVIDENCE_RANK = {name: i for i, (name, _) in enumerate(_EVIDENCE)} | {"none": len(_EVIDENCE)}


def evidence_class(codes: list[str]) -> str:
    for name, members in _EVIDENCE:
        if any(c in members for c in codes):
            return name
    return "none"


# ---------------------------------------------------------------- labels ----

_SUPER_RULES = [
    (re.compile(r"^(?:conotoxin|conopeptide)\s+(\S+)\s+superfamily$", re.I), "conopeptide"),
    (re.compile(r"^(?:conotoxin|conopeptide)\s+(\S+)\s+family$", re.I), "conopeptide"),
    # Bare superfamily codes as UniProt sometimes writes them: "O2 superfamily"
    (re.compile(r"^([A-Z][A-Za-z0-9]{0,3}(?:-like)?)\s+superfamily$"), "conopeptide"),
]
_FAMILY_RULE = re.compile(r"^(.+?)\s+(?:sub)?family$", re.I)


@dataclass
class Label:
    raw_superfamily: str = ""
    raw_family: str = ""
    superfamily: str = ""
    family: str = ""
    scope: str = "unlabelled"   # conopeptide | reference | other | exclude | unlabelled | unmapped
    source: str = ""            # rule | override | ""


def split_similarity(text: str) -> tuple[str, str]:
    """'Belongs to the conotoxin C superfamily. Consomatin family.' ->
    ('conotoxin C superfamily', 'Consomatin family')."""
    text = text.strip()
    m = re.match(r"^Belongs to the (.+?)(?:\.\s+(.+?))?\.?$", text, re.S)
    if not m:
        return text.rstrip("."), ""
    return m.group(1).strip(), (m.group(2) or "").strip().rstrip(".")


@dataclass
class LabelResolver:
    overrides: dict[tuple[str, str], dict[str, str]] = field(default_factory=dict)
    entries: dict[str, dict[str, str]] = field(default_factory=dict)
    taxon_scope: frozenset[str] = frozenset()
    unmapped: dict[tuple[str, str], int] = field(default_factory=dict)

    @classmethod
    def from_config(cls, cfg: LabelsConfig) -> LabelResolver:
        table: dict[tuple[str, str], dict[str, str]] = {}
        if cfg.overrides is not None:
            for row in read_tsv(cfg.overrides):
                raw_sf = (row.get("raw_superfamily") or "").strip()
                if not raw_sf or raw_sf.startswith("#"):
                    continue
                scope = (row.get("scope") or "conopeptide").strip()
                if scope not in {"conopeptide", "other", "exclude"}:
                    raise ValueError(f"{cfg.overrides}: invalid scope {scope!r} for {raw_sf!r}")
                key = (raw_sf.lower(), (row.get("raw_family") or "*").strip().lower() or "*")
                table[key] = {
                    "superfamily": (row.get("superfamily") or "").strip(),
                    "family": (row.get("family") or "").strip(),
                    "scope": scope,
                }
        per_entry: dict[str, dict[str, str]] = {}
        if cfg.entry_overrides is not None:
            for row in read_tsv(cfg.entry_overrides):
                acc = (row.get("accession") or "").strip()
                if not acc or acc.startswith("#"):
                    continue
                scope = (row.get("scope") or "conopeptide").strip()
                if scope not in {"conopeptide", "other", "exclude"}:
                    raise ValueError(f"{cfg.entry_overrides}: invalid scope {scope!r} for {acc}")
                per_entry[acc] = {"superfamily": (row.get("superfamily") or "").strip(),
                                  "family": (row.get("family") or "").strip(), "scope": scope}
        return cls(overrides=table, entries=per_entry, taxon_scope=frozenset(cfg.taxon_scope))

    def apply_taxon(self, taxon_family: str, lab: Label) -> tuple[Label, str | None]:
        """Demote conopeptide labels outside the taxonomic scope to 'reference'."""
        if not self.taxon_scope or lab.scope != "conopeptide":
            return lab, None
        if not taxon_family:
            return lab, "no taxonomic family in lineage; taxon scope not applied"
        if taxon_family not in self.taxon_scope:
            lab.scope = "reference"
        return lab, None

    def apply_entry(self, acc: str, lab: Label) -> Label:
        o = self.entries.get(acc)
        if o is None:
            return lab
        if lab.scope == "unmapped":
            n = self.unmapped.get((lab.raw_superfamily, lab.raw_family), 0) - 1
            if n > 0:
                self.unmapped[(lab.raw_superfamily, lab.raw_family)] = n
            else:
                self.unmapped.pop((lab.raw_superfamily, lab.raw_family), None)
        lab.superfamily, lab.family, lab.scope, lab.source = o["superfamily"], o["family"], o["scope"], "entry_override"
        return lab

    def resolve(self, raw_sf: str, raw_fam: str) -> Label:
        lab = Label(raw_superfamily=raw_sf, raw_family=raw_fam)
        if not raw_sf:
            return lab
        for key in ((raw_sf.lower(), raw_fam.lower()), (raw_sf.lower(), "*")):
            if key in self.overrides:
                o = self.overrides[key]
                lab.superfamily, lab.scope, lab.source = o["superfamily"], o["scope"], "override"
                lab.family = o["family"] or _family_name(raw_fam)
                return lab
        for rx, scope in _SUPER_RULES:
            m = rx.match(raw_sf)
            if m:
                lab.superfamily, lab.scope, lab.source = m.group(1), scope, "rule"
                lab.family = _family_name(raw_fam)
                return lab
        # Not a conopeptide pattern and not curated: keep visible, flag it.
        lab.superfamily = _family_name(raw_sf) or raw_sf
        lab.family = _family_name(raw_fam)
        lab.scope = "unmapped"
        self.unmapped[(raw_sf, raw_fam)] = self.unmapped.get((raw_sf, raw_fam), 0) + 1
        return lab


def _family_name(raw: str) -> str:
    if not raw:
        return ""
    m = _FAMILY_RULE.match(raw)
    return (m.group(1) if m else raw).strip()


# --------------------------------------------------------------- regions ----

@dataclass
class Feature:
    type: str
    start: int | None
    end: int | None
    exact: bool
    codes: list[str]
    description: str


def _features(entry: dict, wanted: set[str]) -> list[Feature]:
    out = []
    for f in entry.get("features", []):
        if f.get("type") not in wanted:
            continue
        loc = f.get("location", {})
        s, e = loc.get("start", {}), loc.get("end", {})
        out.append(
            Feature(
                type=f["type"],
                start=s.get("value"),
                end=e.get("value"),
                exact=s.get("modifier") == "EXACT" and e.get("modifier") == "EXACT",
                codes=[ev.get("evidenceCode", "") for ev in f.get("evidences", [])],
                description=f.get("description", ""),
            )
        )
    return sorted(out, key=lambda x: (x.start or 0, x.end or 0))


def _overlaps(a: Feature, b: Feature) -> bool:
    return a.start <= b.end and b.start <= a.end  # type: ignore[operator]


def extract_regions(entry: dict, cfg: RegionsConfig) -> tuple[list[dict], list[str]]:
    """Return (region rows, warnings) for one entry."""
    acc = entry["primaryAccession"]
    seq = entry.get("sequence", {}).get("value", "")
    warnings: list[str] = []
    feats = _features(entry, {"Signal", "Propeptide", "Peptide", "Chain"})

    usable = []
    for f in feats:
        if f.start is None or f.end is None:
            warnings.append(f"{acc}: {f.type} with unknown boundary skipped")
            continue
        if not f.exact and not cfg.allow_uncertain:
            warnings.append(f"{acc}: {f.type} {f.start}-{f.end} has uncertain boundary, skipped")
            continue
        if not (1 <= f.start <= f.end <= len(seq)):
            warnings.append(f"{acc}: {f.type} {f.start}-{f.end} outside sequence, skipped")
            continue
        usable.append(f)

    sig = [f for f in usable if f.type == "Signal"]
    props = [f for f in usable if f.type == "Propeptide"]
    mats = [f for f in usable if f.type == "Peptide"]
    mat_source = "Peptide"
    if not mats:
        blockers = sig + props
        mats = [f for f in usable if f.type == "Chain" and not any(_overlaps(f, b) for b in blockers)]
        mat_source = "Chain"

    assigned: list[tuple[str, Feature, str]] = [("SIG", f, "Signal") for f in sig[:1]]
    if mats:
        first, last = mats[0].start, mats[-1].end
        before = [p for p in props if p.end < first]  # type: ignore[operator]
        for p in props:
            if any(p is b for b in before):
                region = "PRO" if p is before[-1] else "INTER"
            elif p.start > last:  # type: ignore[operator]
                region = "POST"
            else:
                region = "INTER"
            assigned.append((region, p, "Propeptide"))
        assigned += [("MAT", m, mat_source) for m in mats]
    else:
        assigned += [("INTER", p, "Propeptide") for p in props]

    rows = []
    counters: dict[str, int] = {}
    for region, f, source in sorted(assigned, key=lambda t: t[1].start):  # type: ignore[arg-type,return-value]
        sub = seq[f.start - 1 : f.end]  # type: ignore[operator]
        if len(sub) < cfg.min_length:
            warnings.append(f"{acc}: {region} {f.start}-{f.end} shorter than {cfg.min_length} aa, skipped")
            continue
        counters[region] = counters.get(region, 0) + 1
        rows.append(
            {
                "accession": acc,
                "region": region,
                "index": counters[region],
                "region_id": f"{acc}|{region}|{counters[region]}",
                "start": f.start,
                "end": f.end,
                "length": len(sub),
                "source_feature": source,
                "exact": f.exact,
                "evidence": evidence_class(f.codes),
                "evidence_codes": sorted(set(c for c in f.codes if c)),
                "framework": cysteine_framework(sub) if region == "MAT" else "",
                "sequence": sub,
            }
        )
    return rows, warnings


# ---------------------------------------------------------------- entries ---

def _similarity_text(entry: dict) -> str:
    for c in entry.get("comments", []):
        if c.get("commentType") == "SIMILARITY":
            texts = c.get("texts") or []
            if texts:
                return texts[0].get("value", "")
    return ""


def _protein_name(entry: dict) -> str:
    pd = entry.get("proteinDescription", {})
    for key in ("recommendedName", "submissionNames"):
        v = pd.get(key)
        if isinstance(v, list):
            v = v[0] if v else None
        if v:
            return v.get("fullName", {}).get("value", "")
    return ""


def _taxon_family(lineage: list) -> str:
    names = [x if isinstance(x, str) else x.get("scientificName", "") for x in lineage]
    fams = [n for n in names if n.endswith("idae")]
    return fams[-1] if fams else ""


def parse_entry(entry: dict, resolver: LabelResolver, rcfg: RegionsConfig) -> tuple[dict, list[dict], list[str]]:
    acc = entry["primaryAccession"]
    org = entry.get("organism", {})
    audit = entry.get("entryAudit", {})
    seq = entry.get("sequence", {}).get("value", "")
    flag = entry.get("proteinDescription", {}).get("flag", "") or ""
    raw_sf, raw_fam = split_similarity(_similarity_text(entry))
    taxon_family = _taxon_family(org.get("lineage", []))
    label, taxon_warning = resolver.apply_taxon(taxon_family, resolver.resolve(raw_sf, raw_fam))
    label = resolver.apply_entry(acc, label)
    regions, warnings = extract_regions(entry, rcfg)
    if taxon_warning:
        warnings.append(f"{acc}: {taxon_warning}")
    present = {r["region"] for r in regions}
    row = {
        "accession": acc,
        "entry_name": entry.get("uniProtkbId", ""),
        "protein_name": _protein_name(entry),
        "organism": org.get("scientificName", ""),
        "taxon_id": org.get("taxonId", ""),
        "taxon_family": taxon_family,
        "raw_superfamily": label.raw_superfamily,
        "raw_family": label.raw_family,
        "superfamily": label.superfamily,
        "family": label.family,
        "scope": label.scope,
        "label_source": label.source,
        "length": len(seq),
        "fragment": "fragment" in flag.lower(),
        "complete_precursor": {"SIG", "MAT"} <= present and seq.startswith("M") and "fragment" not in flag.lower(),
        "regions": sorted(present, key=REGIONS.index),
        "first_public": audit.get("firstPublicDate", ""),
        "last_update": audit.get("lastAnnotationUpdateDate", ""),
        "entry_version": audit.get("entryVersion", ""),
        "sequence_version": audit.get("sequenceVersion", ""),
        "sequence": seq,
    }
    return row, regions, warnings


ENTRY_FIELDS = [
    "accession", "entry_name", "protein_name", "organism", "taxon_id", "taxon_family",
    "superfamily", "family", "scope", "label_source", "raw_superfamily", "raw_family",
    "length", "fragment", "complete_precursor", "regions", "first_public", "last_update",
    "entry_version", "sequence_version", "sequence",
]
REGION_FIELDS = [
    "region_id", "accession", "region", "index", "start", "end", "length", "source_feature",
    "exact", "evidence", "evidence_codes", "framework", "sequence",
]


def unmapped_rows(resolver: LabelResolver) -> list[dict]:
    return [
        {"raw_superfamily": sf, "raw_family": fam, "entries": n}
        for (sf, fam), n in sorted(resolver.unmapped.items(), key=lambda kv: -kv[1])
    ]
