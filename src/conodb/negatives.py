"""Negatives: sequences that are not conopeptides, for calibration and the
"none of the known families" class.

Three sources:
  uniprot      a second UniProt query (by default non-venom Conus proteins and
               secreted gastropod proteins, toxins excluded): real proteins,
               most with signal peptides, as conodictor will see them
  venom_other  venom proteins in the main set with scope 'other' (CRISP,
               PLA2, actinoporin...): hard negatives from the venom duct
  shuffled     conopeptide precursors with residues shuffled within each
               region and cysteine positions kept: tests that profiles
               recognise sequence, not just cysteine patterns

QC: real negatives at least `max_identity` identical to any conopeptide entry
are removed (likely unannotated toxins).

Roles: shuffled negatives take their source precursor's side (train or
test), so nothing derived from a test precursor is ever trained on. Real
negatives are clustered at cluster.group_identity and whole groups go to
test with partition.test_fraction, stable across releases like the
positives (inherited from the previous release, new groups by hash order).
"""

from __future__ import annotations

import hashlib
import logging
import random
from collections import Counter, defaultdict
from pathlib import Path

from .cluster import cluster_levels, max_identity
from .config import Config, SourceConfig
from .io import iter_jsonl_gz, read_json, write_json
from .parse import _taxon_family
from .uniprot import ReleaseInfo, fetch_all

log = logging.getLogger(__name__)

NEGATIVE_FIELDS = [
    "id",
    "source",
    "accession",
    "organism",
    "taxon_family",
    "length",
    "role",
    "group",
    "source_role",
    "nearest_conopeptide",
    "max_identity",
    "reason",
    "sequence",
]
_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _h(seed: str, key: str) -> str:
    return hashlib.sha256(f"{seed}|{key}".encode()).hexdigest()


# ----------------------------------------------------------------- fetch ----


def fetch_negatives(
    cfg: Config, release: ReleaseInfo, raw: Path, offline: bool
) -> dict:
    """Download the negatives query once per (query, UniProt release)."""
    meta_path = raw.with_suffix(".meta.json")
    want = {"query": cfg.negatives.query, "release": release.release}
    if raw.exists() and meta_path.exists() and read_json(meta_path) == want:
        return {"downloaded": False}
    if offline:
        log.warning(
            "Offline build: no UniProt negatives (only venom 'other' and shuffled precursors)"
        )
        return {"downloaded": False, "offline": True}
    src = cfg.source
    q = SourceConfig(
        query=cfg.negatives.query,
        base_url=src.base_url,
        page_size=src.page_size,
        timeout=src.timeout,
        retries=src.retries,
    )
    info = fetch_all(q, raw)
    if info.release != release.release:
        raise RuntimeError(
            f"UniProt release changed during the build ({release.release} -> {info.release}); rerun"
        )
    write_json(meta_path, want)
    return {"downloaded": True, "entries": info.total}


# --------------------------------------------------------------- shuffle ----


def shuffle_precursor(
    seq: str, regions: list[tuple[int, int]], rng: random.Random
) -> str:
    """Shuffle residues within each region (and within the stretches between
    regions), keeping every cysteine at its position."""
    seq = seq.upper()
    bounds = sorted((s - 1, e) for s, e in regions)
    cuts = sorted(
        {0, len(seq)} | {b for s, e in bounds for b in (s, e) if 0 <= b <= len(seq)}
    )
    out = list(seq)
    for a, b in zip(cuts, cuts[1:], strict=False):
        idx = [i for i in range(a, b) if seq[i] != "C"]
        vals = [seq[i] for i in idx]
        rng.shuffle(vals)
        for i, v in zip(idx, vals, strict=True):
            out[i] = v
    return "".join(out)


# ----------------------------------------------------------------- build ----


def build_negatives(
    cfg: Config,
    entries: list[dict],
    regions: list[dict],
    partition: list[dict],
    raw: Path | None,
    previous: list[dict] | None,
) -> tuple[list[dict], dict]:
    nc = cfg.negatives
    pc = cfg.partition
    part = {r["accession"]: r for r in partition}
    main_acc = {e["accession"] for e in entries}
    rows: list[dict] = []

    # real negatives
    if raw is not None and raw.exists():
        for e in iter_jsonl_gz(raw):
            acc = e["primaryAccession"]
            seq = e.get("sequence", {}).get("value", "").upper()
            org = e.get("organism", {})
            row = {
                "id": acc,
                "source": "uniprot",
                "accession": acc,
                "organism": org.get("scientificName", ""),
                "taxon_family": _taxon_family(org.get("lineage", [])),
                "length": len(seq),
                "sequence": seq,
            }
            if acc in main_acc:
                row |= {"role": "removed", "reason": "also in the conopeptide query"}
            elif not nc.min_length <= len(seq) <= nc.max_length:
                row |= {
                    "role": "removed",
                    "reason": f"length outside {nc.min_length}-{nc.max_length}",
                }
            elif set(seq) - _AA - {"X", "U", "B", "Z"}:
                row |= {"role": "removed", "reason": "non-protein characters"}
            rows.append(row)
    if nc.include_venom_other:
        for e in entries:
            if e["scope"] == "other":
                rows.append(
                    {
                        "id": e["accession"],
                        "source": "venom_other",
                        "accession": e["accession"],
                        "organism": e["organism"],
                        "taxon_family": e["taxon_family"],
                        "length": int(e["length"]),
                        "sequence": e["sequence"].upper(),
                    }
                )

    # cap on real UniProt negatives: stable choice by accession hash
    if nc.max_entries > 0:
        cand = sorted(
            (r for r in rows if r["source"] == "uniprot" and not r.get("role")),
            key=lambda r: _h(pc.seed, "neg|" + r["id"]),
        )
        for r in cand[nc.max_entries :]:
            r |= {
                "role": "removed",
                "reason": f"over negatives.max_entries ({nc.max_entries})",
            }

    # QC against conopeptides
    cono = [
        (e["accession"], e["sequence"])
        for e in entries
        if e["scope"] == "conopeptide" and e["sequence"]
    ]
    real = [r for r in rows if not r.get("role")]
    nearest = max_identity(
        ((r["id"], r["sequence"]) for r in real), cono, cfg.cluster.threads
    )
    for r in real:
        tid, ident = nearest.get(r["id"], ("", 0.0))
        r["nearest_conopeptide"], r["max_identity"] = tid, round(ident, 4)
        if ident >= nc.max_identity:
            r |= {
                "role": "removed",
                "reason": f">= {nc.max_identity:g} identical to conopeptide {tid}",
            }

    # real negatives: groups and a stable test split
    kept = [r for r in rows if not r.get("role")]
    levels = cluster_levels(
        [(r["id"], r["sequence"]) for r in kept],
        [cfg.cluster.group_identity],
        cfg.cluster.threads,
    )
    group = levels.get(cfg.cluster.group_identity, {})
    prev = {
        r["id"]: r["role"] for r in previous or [] if r.get("role") in ("train", "test")
    }
    by_g: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in kept:
        r["group"] = group[r["id"]]
        by_g[(r["source"], r["group"])].append(r)
    by_source: dict[str, list[tuple[str, list[dict]]]] = defaultdict(list)
    for (src, g), members in by_g.items():
        by_source[src].append((g, members))
    for _src, gs in by_source.items():
        gs.sort(key=lambda gm: _h(pc.seed, min(m["id"] for m in gm[1])))
        role = {}
        for g, members in gs:
            old = {prev[m["id"]] for m in members if m["id"] in prev}
            if "test" in old:
                role[g] = "test"
            elif "train" in old:
                role[g] = "train"
        quota = round(pc.test_fraction * len(gs))
        n_test = sum(v == "test" for v in role.values())
        for g, _ in gs:
            if g not in role:
                role[g] = "test" if n_test < quota else "train"
                n_test += role[g] == "test"
        for g, members in gs:
            for m in members:
                m["role"] = role[g]
                m["reason"] = "inherited" if m["id"] in prev else ""

    # shuffled precursors: one or more per complete precursor, same side as the source
    if nc.shuffles > 0:
        regs = defaultdict(list)
        for r in regions:
            regs[r["accession"]].append((int(r["start"]), int(r["end"])))
        for e in entries:
            p = part.get(e["accession"])
            if (
                not p
                or p["role"] not in ("train", "test", "temporal_test")
                or e["complete_precursor"] != "true"
            ):
                continue
            side = "train" if p["role"] == "train" else "test"
            for k in range(nc.shuffles):
                rng = random.Random(_h(pc.seed, f"shuffle|{e['accession']}|{k}"))
                seq = shuffle_precursor(e["sequence"], regs[e["accession"]], rng)
                rows.append(
                    {
                        "id": f"SHUF-{e['accession']}-{k + 1}",
                        "source": "shuffled",
                        "accession": e["accession"],
                        "organism": e["organism"],
                        "taxon_family": e["taxon_family"],
                        "length": len(seq),
                        "role": side,
                        "group": f"SHUF-{p['group']}",
                        "source_role": p["role"],
                        "reason": "",
                        "sequence": seq,
                    }
                )

    rows.sort(key=lambda r: (r["source"], r["id"]))
    stats = {
        "by_source_role": {
            f"{s}:{ro}": n
            for (s, ro), n in sorted(
                Counter((r["source"], r["role"]) for r in rows).items()
            )
        },
        "removed_similar_to_conopeptide": sum(
            r["role"] == "removed" and r["reason"].startswith(">=") for r in rows
        ),
    }
    return rows, stats
