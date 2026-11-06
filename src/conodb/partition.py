"""Partition: held-out test set and cross-validation folds, stable across releases.

Units are entry groups (full-sequence clusters at cluster.group_identity):
a group is never split between training and test.

Per superfamily, after removing excluded entries and recomputing tiers:
  tier 1  ~test_fraction of the groups -> 'test'; the rest -> 'train' with
          grouped folds 1..k (balanced by number of groups)
  tier 2  'train' with fold 'logo' (leave-one-group-out cross-validation)
  tier 3  'reference' (reference sequences only; not a class, not evaluated)
Optionally, groups first made public on/after partition.temporal_cutoff (in
tier 1-2 superfamilies) -> 'temporal_test'.

Stability: when a previous release exists, every group inherits the role and
fold its members had (a group containing any former test entry stays test).
Only new groups are assigned, by a hash of their smallest accession, so a
database update adds data without reshuffling the benchmark. Groups are moved
from train to test only when a superfamily's test quota cannot be met with
new groups (reported in partition.md).
"""

from __future__ import annotations

import hashlib
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .audit import tier_of
from .config import Config

log = logging.getLogger(__name__)

PARTITION_FIELDS = ["accession", "superfamily", "family", "scope", "tier", "group", "role", "fold", "inherited", "reason"]
SUMMARY_FIELDS = ["superfamily", "tier", "groups", "entries", "train_groups", "test_groups", "temporal_test_groups",
                  "train_entries", "test_entries", "temporal_test_entries", "fold_groups"]

_SCOPE_ROLE = {"other": "negative", "reference": "reference", "unlabelled": "unlabelled",
               "unmapped": "unmapped", "exclude": "excluded"}


def _h(seed: str, key: str) -> str:
    return hashlib.sha256(f"{seed}|{key}".encode()).hexdigest()


@dataclass
class Group:
    gid: str
    members: list[dict]
    superfamily: str = ""
    mixed: bool = False
    key: str = ""                 # smallest accession: stable identity across releases
    first_public: str = ""
    role: str = ""
    fold: str = ""
    inherited: bool = False
    reason: str = ""
    notes: list[str] = field(default_factory=list)


@dataclass
class PartitionResult:
    rows: list[dict]
    summary: list[dict]
    excluded: list[dict]
    reassigned: list[str]
    tier_changes: list[dict]
    mixed_groups: list[str]


def run_partition(cfg: Config, entries: list[dict], groups_tsv: list[dict], conflicts: list[dict],
                  audit_tiers: dict[str, str], previous: list[dict] | None) -> PartitionResult:
    pc = cfg.partition
    group_of = {r["accession"]: r["group"] for r in groups_tsv}

    # 1. exclusions: accessions under review in label conflicts
    excluded: dict[str, str] = {}
    if pc.exclude_conflicts:
        for c in conflicts:
            for acc in filter(None, c.get("review_accessions", "").split(";")):
                excluded.setdefault(acc, f"label conflict ({c['severity']}: {c['labels']})")

    # 2. eligible entries -> groups
    eligible = [e for e in entries if e["scope"] == "conopeptide" and e["accession"] not in excluded]
    by_gid: dict[str, list[dict]] = defaultdict(list)
    for e in eligible:
        by_gid[group_of[e["accession"]]].append(e)
    groups: list[Group] = []
    for gid, members in sorted(by_gid.items()):
        labels = Counter(m["superfamily"] for m in members)
        top = max(labels.values())
        sf = sorted(k for k, v in labels.items() if v == top)[0]
        dates = [m["first_public"] for m in members if m["first_public"]]
        groups.append(Group(gid=gid, members=members, superfamily=sf, mixed=len(labels) > 1,
                            key=min(m["accession"] for m in members), first_public=min(dates) if dates else ""))

    # 3. tiers after exclusions
    n_groups = Counter(g.superfamily for g in groups)
    tiers = {sf: tier_of(n, cfg) for sf, n in n_groups.items()}
    tier_changes = [{"superfamily": sf, "audit_tier": audit_tiers.get(sf, ""), "partition_tier": t}
                    for sf, t in sorted(tiers.items()) if str(audit_tiers.get(sf, "")) != str(t)]

    # 4. inheritance from the previous release
    prev = {r["accession"]: r for r in previous or []}
    prev_k = {r["fold"] for r in previous or [] if r.get("fold", "").isdigit()}
    folds_ok = bool(prev_k) and max(int(f) for f in prev_k) <= pc.folds
    if previous and prev_k and not folds_ok:
        log.warning("partition.folds changed; training folds are reassigned")
    for g in groups:
        roles = {prev[m["accession"]]["role"] for m in g.members if m["accession"] in prev}
        if "test" in roles:
            g.role, g.inherited = "test", True
        elif "temporal_test" in roles:
            g.role, g.inherited = "temporal_test", True
        elif "train" in roles:
            g.role, g.inherited = "train", True
            folds = Counter(prev[m["accession"]]["fold"] for m in g.members
                            if m["accession"] in prev and prev[m["accession"]]["fold"].isdigit())
            if folds and folds_ok:
                g.fold = sorted(folds.items(), key=lambda kv: (-kv[1], int(kv[0])))[0][0]

    # 5. roles by tier
    reassigned: list[str] = []
    by_sf: dict[str, list[Group]] = defaultdict(list)
    for g in groups:
        by_sf[g.superfamily].append(g)
    for sf, gs in sorted(by_sf.items()):
        tier = tiers[sf]
        gs.sort(key=lambda g: _h(pc.seed, g.key))
        if tier == 3:
            for g in gs:
                g.role, g.fold, g.reason = "reference", "", "tier 3"
            continue
        if pc.temporal_cutoff:
            for g in gs:
                if not g.inherited and g.first_public and g.first_public >= pc.temporal_cutoff:
                    g.role, g.reason = "temporal_test", f"first public >= {pc.temporal_cutoff}"
        if tier == 2:
            for g in gs:
                if g.role != "temporal_test":
                    g.role, g.fold = "train", "logo"
            continue
        # tier 1: test quota over non-temporal groups
        pool = [g for g in gs if g.role != "temporal_test"]
        quota = max(1, round(pc.test_fraction * len(pool)))
        n_test = sum(g.role == "test" for g in pool)
        for g in pool:                       # new groups first, in hash order
            if n_test >= quota:
                break
            if not g.role:
                g.role, g.reason = "test", "new group"
                n_test += 1
        for g in pool:                       # then, only if unavoidable, former training groups
            if n_test >= quota:
                break
            if g.role == "train":
                g.role, g.fold, g.reason = "test", "", "moved from train to meet test quota"
                reassigned.append(f"{sf}: group {g.key} ({len(g.members)} entries)")
                n_test += 1
        # training folds, balanced by group count
        train = [g for g in pool if g.role in ("", "train")]
        sizes = Counter(g.fold for g in train if g.fold)
        for g in train:
            g.role = "train"
            if not g.fold:
                g.fold = str(min(range(1, pc.folds + 1), key=lambda f: (sizes[str(f)], f)))
                sizes[g.fold] += 1

    # 6. rows for every entry
    rows = []
    gmap = {m["accession"]: g for g in groups for m in g.members}
    for e in entries:
        acc = e["accession"]
        base = {"accession": acc, "superfamily": e["superfamily"], "family": e["family"], "scope": e["scope"],
                "group": group_of.get(acc, "")}
        if acc in gmap:
            g = gmap[acc]
            rows.append(base | {"tier": tiers[g.superfamily], "role": g.role, "fold": g.fold, "inherited": g.inherited,
                                "reason": g.reason + (f"; group labelled {g.superfamily}" if g.mixed else "")})
        elif acc in excluded:
            rows.append(base | {"tier": "", "role": "excluded", "fold": "", "inherited": False, "reason": excluded[acc]})
        else:
            rows.append(base | {"tier": "", "role": _SCOPE_ROLE.get(e["scope"], e["scope"]), "fold": "",
                                "inherited": False, "reason": f"scope {e['scope']}"})
    rows.sort(key=lambda r: r["accession"])

    # 7. integrity: a group never spans train and a test role
    split = defaultdict(set)
    for r in rows:
        if r["role"] in ("train", "test", "temporal_test"):
            split[r["group"]].add("train" if r["role"] == "train" else "test")
    leaks = [g for g, s in split.items() if len(s) > 1]
    if leaks:
        raise RuntimeError(f"Partition leak: groups in both train and test: {leaks[:5]}")

    # 8. summary
    summary = []
    for sf, gs in sorted(by_sf.items(), key=lambda kv: (tiers[kv[0]], -len(kv[1]), kv[0])):
        def cnt(role, what, gs=gs):
            sel = [g for g in gs if g.role == role]
            return len(sel) if what == "g" else sum(len(g.members) for g in sel)
        fold_groups = Counter(g.fold for g in gs if g.role == "train" and g.fold.isdigit())
        summary.append({
            "superfamily": sf, "tier": tiers[sf], "groups": len(gs), "entries": sum(len(g.members) for g in gs),
            "train_groups": cnt("train", "g"), "test_groups": cnt("test", "g"), "temporal_test_groups": cnt("temporal_test", "g"),
            "train_entries": cnt("train", "e"), "test_entries": cnt("test", "e"), "temporal_test_entries": cnt("temporal_test", "e"),
            "fold_groups": "/".join(str(fold_groups[str(f)]) for f in range(1, pc.folds + 1)) if tiers[sf] == 1
            else ("logo" if tiers[sf] == 2 else ""),
        })
    excl_rows = [{"accession": a, "reason": r} for a, r in sorted(excluded.items())]
    mixed = [f"{g.key} ({', '.join(sorted({m['superfamily'] for m in g.members}))})" for g in groups if g.mixed]
    return PartitionResult(rows, summary, excl_rows, reassigned, tier_changes, mixed)


def render_partition(cfg: Config, release: str, res: PartitionResult, previous: str | None) -> str:
    from .audit import _table

    pc = cfg.partition
    roles = Counter(r["role"] for r in res.rows)
    inherited = sum(r["inherited"] is True for r in res.rows)
    lines = [
        f"# Partition: UniProt {release}",
        "",
        f"- Units: entry groups at {cfg.cluster.group_identity:g} identity; a group is never split.",
        f"- Tier 1: {pc.test_fraction:.0%} of groups held out as test, {pc.folds} grouped folds on the rest. "
        "Tier 2: leave-one-group-out. Tier 3: reference only.",
        f"- Temporal test: {'groups first public on/after ' + pc.temporal_cutoff if pc.temporal_cutoff else 'off'}.",
        f"- Seed: `{pc.seed}`. Previous release: {previous or 'none (fresh partition)'}"
        + (f"; {inherited} entries kept their previous assignment." if previous else "."),
        "",
        "Entries by role: " + ", ".join(f"{k} {v}" for k, v in sorted(roles.items(), key=lambda kv: -kv[1])),
        "",
        "## Per superfamily",
        "",
        "`fold_groups`: training groups in each fold (tier 1).",
        "",
        _table(res.summary, ["superfamily", "tier", "groups", "entries", "train_groups", "test_groups",
                             "temporal_test_groups", "train_entries", "test_entries", "fold_groups"]),
        f"## Excluded entries ({len(res.excluded)})",
        "",
        _table(res.excluded, ["accession", "reason"]),
        "## Tier changes after exclusions",
        "",
        _table(res.tier_changes, ["superfamily", "audit_tier", "partition_tier"]),
        "## Groups mixing superfamily labels",
        "",
        ("\n".join(f"- {m}" for m in res.mixed_groups) + "\n") if res.mixed_groups else "_none_\n",
        "## Groups moved from train to test",
        "",
        ("\n".join(f"- {m}" for m in res.reassigned) + "\n") if res.reassigned else "_none_\n",
    ]
    return "\n".join(lines) + "\n"
