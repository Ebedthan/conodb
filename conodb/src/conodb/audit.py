"""Audit: how much independent evidence each (super)family has, per region.

Produces the tables that decide support tiers and the partition design,
flags label inconsistencies, and compares with the previous release.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .config import Config
from .parse import EVIDENCE_RANK, PROFILE_REGIONS

ALL = "*"


def tier_of(groups: int, cfg: Config) -> int:
    if groups >= cfg.tiers.tier1_min_groups:
        return 1
    if groups >= cfg.tiers.tier2_min_groups:
        return 2
    return 3


def _keys(entry: dict) -> list[tuple[str, str]]:
    sf = entry["superfamily"] or "(unlabelled)"
    keys = [(sf, ALL)]
    if entry["family"]:
        keys.append((sf, entry["family"]))
    return keys


@dataclass
class AuditResult:
    families: list[dict]
    region_support: list[dict]
    conflicts: list[dict]
    unlabelled: dict
    taxa: list[dict]


def run_audit(cfg: Config, entries: list[dict], regions: list[dict], clusters: list[dict], groups: dict[str, str]) -> AuditResult:
    idents = sorted(cfg.cluster.identities, reverse=True)
    by_acc = {e["accession"]: e for e in entries}
    cl: dict[tuple[str, float], str] = {(c["item_id"], float(c["identity"])): c["cluster"] for c in clusters}

    fam_entries: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for e in entries:
        if e["scope"] == "exclude":
            continue
        for sf, fam in _keys(e):
            fam_entries[(sf, fam, e["scope"])].append(e)

    full_cl = {t: {c["item_id"]: c["cluster"] for c in clusters if c["kind"] == "FULL" and float(c["identity"]) == t}
               for t in idents}
    families = []
    for (sf, fam, scope), members in sorted(fam_entries.items()):
        n_groups = len({groups[m["accession"]] for m in members})
        by_t = {f"groups_{t:g}": len({full_cl[t][m["accession"]] for m in members}) for t in idents}
        is_conopeptide = scope == "conopeptide"
        families.append(
            {
                "superfamily": sf,
                "family": fam,
                "scope": scope,
                "entries": len(members),
                "entry_groups": n_groups,
                "species": len({m["organism"] for m in members}),
                "complete_precursors": sum(m["complete_precursor"] == "true" for m in members),
                "fragments": sum(m["fragment"] == "true" for m in members),
                "tier": tier_of(n_groups, cfg) if is_conopeptide else "",
                "first_public_min": min((m["first_public"] for m in members if m["first_public"]), default=""),
                "first_public_max": max((m["first_public"] for m in members if m["first_public"]), default=""),
                **by_t,
            }
        )

    reg_by_key: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for r in regions:
        e = by_acc.get(r["accession"])
        if e is None or e["scope"] != "conopeptide" or r["region"] not in PROFILE_REGIONS:
            continue
        for sf, fam in _keys(e):
            reg_by_key[(sf, fam, r["region"])].append(r)

    support = []
    for (sf, fam, region), rows in sorted(reg_by_key.items(), key=lambda kv: (kv[0][0], kv[0][1], PROFILE_REGIONS.index(kv[0][2]))):
        row = {
            "superfamily": sf,
            "family": fam,
            "region": region,
            "sequences": len(rows),
            "unique": len({r["sequence"] for r in rows}),
            "entries": len({r["accession"] for r in rows}),
            "exact": sum(r["exact"] == "true" for r in rows),
            "experimental": sum(EVIDENCE_RANK[r["evidence"]] == 0 for r in rows),
            "median_length": sorted(int(r["length"]) for r in rows)[len(rows) // 2],
            "frameworks": len({r["framework"] for r in rows if r["framework"]}) if region == "MAT" else "",
        }
        for t in idents:
            row[f"clusters_{t:g}"] = len({cl[(r["region_id"], t)] for r in rows})
        support.append(row)

    # Label QC: near-identical sequences (highest threshold) carrying
    # different conopeptide superfamily labels. Short regions are ignored:
    # they are identical across families by chance (e.g. 'GKR').
    top = idents[0]
    kinds = set(cfg.audit.conflict_kinds)
    min_len = cfg.audit.conflict_min_length
    length = {r["region_id"]: int(r["length"]) for r in regions} | {e["accession"]: int(e["length"]) for e in entries}
    members: dict[tuple[str, str], set[str]] = defaultdict(set)
    for c in clusters:
        if float(c["identity"]) != top or c["kind"] not in kinds or length.get(c["item_id"], 0) < min_len:
            continue
        acc = c["item_id"].split("|")[0]
        e = by_acc.get(acc)
        if e and e["scope"] == "conopeptide":
            members[(c["kind"], c["cluster"])].add(acc)
    # One row per set of conflicting entries, listing every sequence kind in
    # which they are near-identical, with the minority labels to review.
    merged: dict[frozenset[str], dict] = {}
    for (kind, rep), accs in sorted(members.items()):
        labels = defaultdict(list)
        for a in sorted(accs):
            labels[by_acc[a]["superfamily"]].append(a)
        if len(labels) <= 1:
            continue
        top_n = max(len(v) for v in labels.values())
        leaders = sorted(k for k, v in labels.items() if len(v) == top_n)
        # A tie has no majority: every tied label is under review.
        majority = leaders[0] if len(leaders) == 1 else ""
        row = merged.setdefault(frozenset(accs), {
            "identity": top,
            "representative": rep.split("|")[0],
            "labels": ", ".join(f"{k} x{len(v)}" for k, v in sorted(labels.items(), key=lambda kv: -len(kv[1]))),
            "majority": majority,
            "review": sorted(f"{a}:{k}" for k, v in labels.items() if k != majority for a in v),
            "review_accessions": sorted(a for k, v in labels.items() if k != majority for a in v),
            "accessions": sorted(f"{a}:{by_acc[a]['superfamily']}" for a in accs),
            "kind": [],
        })
        row["kind"].append(kind)
    order = {k: i for i, k in enumerate(CONFLICT_KIND_ORDER)}
    for row in merged.values():
        row["kind"].sort(key=order.get)
        row["severity"] = row["kind"][0]
    conflicts = sorted(merged.values(), key=lambda r: (order[r["severity"]], r["representative"]))
    taxa = defaultdict(lambda: defaultdict(int))
    for e in entries:
        taxa[e["taxon_family"] or "(unknown)"][e["scope"]] += 1
    taxa_rows = [{"taxon_family": t, "entries": sum(c.values()),
                  "by_scope": ", ".join(f"{k} {v}" for k, v in sorted(c.items(), key=lambda kv: -kv[1]))}
                 for t, c in sorted(taxa.items(), key=lambda kv: -sum(kv[1].values()))]
    return AuditResult(families, support, conflicts, _unlabelled_summary(entries, regions), taxa_rows)


def _unlabelled_summary(entries: list[dict], regions: list[dict]) -> dict:
    """What the entries without a UniProt family look like: are they usable
    precursors (future soft labels / predictions) or bare mature peptides?"""
    unl = [e for e in entries if e["scope"] == "unlabelled"]
    regs = defaultdict(set)
    for r in regions:
        regs[r["accession"]].add(r["region"])
    stems = defaultdict(int)
    for e in unl:
        name = e["protein_name"] or "(no name)"
        stems[name.split()[0].rstrip(",")] += 1
    return {
        "entries": len(unl),
        "complete_precursors": sum(e["complete_precursor"] == "true" for e in unl),
        "with_signal": sum("SIG" in regs[e["accession"]] for e in unl),
        "mature_only": sum(regs[e["accession"]] == {"MAT"} for e in unl),
        "fragments": sum(e["fragment"] == "true" for e in unl),
        "name_stems": sorted(stems.items(), key=lambda kv: (-kv[1], kv[0]))[:12],
    }


CONFLICT_KIND_ORDER = ["FULL", "SIG", "MAT", "PRO", "POST"]


def family_fields(cfg: Config) -> list[str]:
    idents = sorted(cfg.cluster.identities, reverse=True)
    return ["superfamily", "family", "scope", "tier", "entries", "entry_groups", "species",
            "complete_precursors", "fragments", "first_public_min", "first_public_max"] + [f"groups_{t:g}" for t in idents]


def support_fields(cfg: Config) -> list[str]:
    idents = sorted(cfg.cluster.identities, reverse=True)
    return ["superfamily", "family", "region", "sequences", "unique", "entries", "exact", "experimental",
            "median_length", "frameworks"] + [f"clusters_{t:g}" for t in idents]


CONFLICT_FIELDS = ["severity", "kind", "majority", "labels", "review", "review_accessions", "accessions", "identity"]


# -------------------------------------------------------------------- diff --

def diff_entries(old: list[dict], new: list[dict]) -> list[dict]:
    o = {e["accession"]: e for e in old}
    n = {e["accession"]: e for e in new}
    out = []
    for acc in sorted(n.keys() - o.keys()):
        out.append({"accession": acc, "change": "added", "detail": n[acc]["superfamily"]})
    for acc in sorted(o.keys() - n.keys()):
        out.append({"accession": acc, "change": "removed", "detail": o[acc]["superfamily"]})
    for acc in sorted(o.keys() & n.keys()):
        a, b = o[acc], n[acc]
        if a["sequence"] != b["sequence"]:
            out.append({"accession": acc, "change": "sequence", "detail": f"v{a['sequence_version']}->v{b['sequence_version']}"})
        if (a["superfamily"], a["family"], a["scope"]) != (b["superfamily"], b["family"], b["scope"]):
            out.append({"accession": acc, "change": "label",
                        "detail": f"{a['superfamily']}/{a['family'] or '-'} -> {b['superfamily']}/{b['family'] or '-'}"})
        elif a["regions"] != b["regions"]:
            out.append({"accession": acc, "change": "regions", "detail": f"{a['regions']} -> {b['regions']}"})
    return out


def diff_tiers(old: list[dict], new: list[dict]) -> list[dict]:
    o = {(f["superfamily"], f["family"], f.get("scope", "")): f for f in old}
    out = []
    for f in new:
        if f.get("scope") != "conopeptide":
            continue
        prev = o.get((f["superfamily"], f["family"], f.get("scope", "")))
        if prev is None:
            out.append({"superfamily": f["superfamily"], "family": f["family"], "old": "new", "new": f["tier"]})
        elif str(prev["tier"]) != str(f["tier"]):
            out.append({"superfamily": f["superfamily"], "family": f["family"], "old": prev["tier"], "new": f["tier"]})
    return out


DIFF_FIELDS = ["accession", "change", "detail"]


# ------------------------------------------------------------------ report --

def _table(rows: list[dict], cols: list[str]) -> str:
    if not rows:
        return "_none_\n"
    head = "| " + " | ".join(cols) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
    return head + "".join("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |\n" for r in rows)


def render_report(
    cfg: Config,
    release: dict,
    result: AuditResult,
    unmapped: list[dict],
    warn_counts: dict[str, int],
    entry_diff: list[dict] | None,
    tier_diff: list[dict] | None,
    previous: str | None,
) -> str:
    idents = sorted(cfg.cluster.identities, reverse=True)
    sf_rows = [f for f in result.families if f["family"] == ALL]
    cono = sorted((f for f in sf_rows if f["scope"] == "conopeptide"), key=lambda f: (-int(f["entry_groups"]), f["superfamily"]))
    others = [f for f in sf_rows if f["scope"] != "conopeptide"]
    sub = [f for f in result.families if f["family"] != ALL and f["scope"] == "conopeptide"]
    tiers = defaultdict(int)
    for f in cono:
        tiers[f["tier"]] += 1
    t = cfg.tiers

    # Sensitivity of the tiers to the grouping threshold.
    sens = []
    for th in idents:
        counts = defaultdict(int)
        for f in cono:
            counts[tier_of(int(f[f"groups_{th:g}"]), cfg)] += 1
        sens.append({"identity": f"{th:g}" + (" (current)" if th == cfg.cluster.group_identity else ""),
                     "tier 1": counts[1], "tier 2": counts[2], "tier 3": counts[3]})

    u = result.unlabelled
    kinds = defaultdict(int)
    for c in result.conflicts:
        kinds[c["severity"]] += 1

    lines = [
        f"# ConoDictor database audit: UniProt {release['release']}",
        "",
        f"- UniProt release: **{release['release']}** {('(' + release['release_date'] + ')') if release.get('release_date') else ''}",
        f"- Query: `{cfg.source.query}`",
        f"- Entries downloaded: **{release['total']}**",
        f"- Conopeptide superfamilies: **{len(cono)}** "
        f"(tier 1: {tiers[1]}, tier 2: {tiers[2]}, tier 3: {tiers[3]})",
        f"- Tiers count entry groups at {cfg.cluster.group_identity:g} identity: "
        f"tier 1 >= {t.tier1_min_groups}, tier 2 >= {t.tier2_min_groups}, tier 3 below.",
        "",
        "## Entries by taxonomic family",
        "",
        _table(result.taxa, ["taxon_family", "entries", "by_scope"]),
        "## Conopeptide superfamilies",
        "",
        "`groups_x` = independent entry groups if full sequences are grouped at identity x.",
        "",
        _table(cono, ["superfamily", "tier", "entries", "entry_groups", "species", "complete_precursors", "fragments"]
               + [f"groups_{th:g}" for th in idents]),
        "### Tiers under each grouping threshold",
        "",
        _table(sens, ["identity", "tier 1", "tier 2", "tier 3"]),
        "## Families within superfamilies",
        "",
        _table(sub, ["superfamily", "family", "tier", "entries", "entry_groups", "species"]),
        "## Other entries (not classifier classes)",
        "",
        "`reference`: conopeptide families outside `labels.taxon_scope`; `other`: venom proteins kept as "
        "future negatives; `unmapped`: awaiting curation.",
        "",
        _table(others, ["superfamily", "scope", "entries", "entry_groups"]),
        "## Entries without a UniProt family",
        "",
        f"{u['entries']} entries: {u['complete_precursors']} complete precursors, {u['with_signal']} with a signal "
        f"peptide, {u['mature_only']} mature peptide only, {u['fragments']} flagged as fragments.",
        "",
        "Most frequent first word of the protein name: "
        + (", ".join(f"{k} ({v})" for k, v in u["name_stems"]) or "-"),
        "",
        "## Labels needing curation",
        "",
        "Raw UniProt family text that no rule or override covers. Add rows to the labels file "
        "with scope `conopeptide`, `other` or `exclude`.",
        "",
        _table(unmapped, ["raw_superfamily", "raw_family", "entries"]),
        f"## Label conflicts ({len(result.conflicts)})",
        "",
        f"Sequences >= {idents[0]:g} identical (regions >= {cfg.audit.conflict_min_length} aa, kinds: "
        f"{', '.join(cfg.audit.conflict_kinds)}) but labelled with different superfamilies, most severe first. "
        "FULL = whole sequences; SIG = signal peptides (superfamilies are defined by signal similarity). "
        "`review` lists the minority labels. By severity: "
        + (", ".join(f"{k} {kinds[k]}" for k in CONFLICT_KIND_ORDER if kinds[k]) or "none") + ".",
        "",
        _table([{**c, "kind": ", ".join(c["kind"]), "review": ", ".join(c["review"])} for c in result.conflicts[:30]],
               ["severity", "kind", "labels", "review"]),
        "## Parsing warnings",
        "",
        _table([{"warning": k, "count": v} for k, v in sorted(warn_counts.items())], ["warning", "count"]),
    ]
    if previous is not None and entry_diff is not None:
        counts = defaultdict(int)
        for d in entry_diff:
            counts[d["change"]] += 1
        lines += [
            f"## Changes since {previous}",
            "",
            ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())) or "No entry changes.",
            "",
            "### Tier changes",
            "",
            _table(tier_diff or [], ["superfamily", "family", "old", "new"]),
        ]
    return "\n".join(lines) + "\n"
