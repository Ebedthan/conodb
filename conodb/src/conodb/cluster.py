"""Greedy, CD-HIT-style identity clustering with pyopal (no external binary).

identity(a, b) = identical aligned residues / length of the shorter sequence,
from a local (Smith-Waterman) alignment. A mature-peptide fragment contained
in a full precursor therefore has identity 1.0 to it, which is what we want
for leakage-aware grouping.

Thresholds are applied hierarchically (highest first, then clustering the
representatives), so clusters are nested across thresholds and the cost
stays low as the database grows.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterable
from multiprocessing.pool import ThreadPool

import pyopal

log = logging.getLogger(__name__)


def _identity(result, qlen: int) -> float:
    return result.alignment.count("M") / max(1, min(qlen, result.target_length))


def greedy_cluster(items: Iterable[tuple[str, str]], identity: float, threads: int = 0, top_k: int = 8) -> dict[str, str]:
    """Map each item id to the id of its cluster representative."""
    uniq: dict[str, list[str]] = {}
    for name, seq in items:
        uniq.setdefault(seq.upper(), []).append(name)
    order = sorted(uniq, key=lambda s: (-len(s), min(uniq[s])))

    reps = pyopal.Database()        # grows in place: no re-encoding per query
    rep_names: list[str] = []
    assign: dict[str, str] = {}
    with ThreadPool(threads or os.cpu_count() or 1) as pool:
        for seq in order:
            names = sorted(uniq[seq])
            best: tuple[int, float] | None = None
            if rep_names:
                # Cheap score pass over all representatives, full alignment
                # only for the top-k candidates.
                # Deterministic: ties broken by representative order, so the
                # same input always gives the same clusters.
                scores = pyopal.align(seq, reps, mode="score", pool=pool, ordered=True)
                cand = [r.target_index for r in sorted(scores, key=lambda r: (-r.score, r.target_index))[:top_k]]
                cand.sort()
                for r in pyopal.align(seq, reps.extract(cand), mode="full", pool=pool, ordered=True):
                    ident = _identity(r, len(seq))
                    if ident >= identity and (best is None or ident > best[1]):
                        best = (cand[r.target_index], ident)
            if best is None:
                reps.append(seq)
                rep_names.append(names[0])
                rep = names[0]
            else:
                rep = rep_names[best[0]]
            for n in names:
                assign[n] = rep
    return assign


def cluster_levels(items: list[tuple[str, str]], identities: list[float], threads: int = 0) -> dict[float, dict[str, str]]:
    """Nested clustering at several thresholds: {threshold: {id: representative}}."""
    seqs = dict(items)
    out: dict[float, dict[str, str]] = {}
    current = {name: name for name in seqs}  # id -> rep at previous level
    for t in sorted(set(identities), reverse=True):
        reps = sorted(set(current.values()))
        rep_assign = greedy_cluster(((r, seqs[r]) for r in reps), t, threads)
        current = {name: rep_assign[rep] for name, rep in current.items()}
        out[t] = dict(current)
        log.debug("  identity %.2f: %d items -> %d clusters", t, len(seqs), len(set(current.values())))
    return out


def max_identity(queries: Iterable[tuple[str, str]], targets: Iterable[tuple[str, str]], threads: int = 0,
                 top_k: int = 8) -> dict[str, tuple[str, float]]:
    """For each query, the most identical target: {query_id: (target_id, identity)}.

    Same identity definition as the clustering (identical residues / length
    of the shorter sequence). Used for negative-set QC and, at evaluation,
    for reporting accuracy by identity to the training set.
    """
    targets = sorted(targets)
    if not targets:
        return {}
    db = pyopal.Database([s.upper() for _, s in targets])
    out: dict[str, tuple[str, float]] = {}
    with ThreadPool(threads or os.cpu_count() or 1) as pool:
        for qid, seq in sorted(queries):
            seq = seq.upper()
            scores = pyopal.align(seq, db, mode="score", pool=pool, ordered=True)
            cand = sorted(r.target_index for r in sorted(scores, key=lambda r: (-r.score, r.target_index))[:top_k])
            best = ("", 0.0)
            for r in pyopal.align(seq, db.extract(cand), mode="full", pool=pool, ordered=True):
                ident = _identity(r, len(seq))
                if ident > best[1]:
                    best = (targets[cand[r.target_index]][0], ident)
            out[qid] = best
    return out
