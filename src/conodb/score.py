"""Score: out-of-fold profile scores for every positive and negative sequence.

Scored sequences: conopeptide entries in train/test/temporal_test, and
negatives in train/test. Every sequence is scored against the final
profiles of every family. A training entry is additionally scored against
its own family's holdout profiles for its unit (its fold, or its group in
tier 2), and in the feature table that holdout score replaces the final one.
Every feature is therefore computed with profiles that never saw the
sequence (test entries and negatives are in no profile at all).

HMM: pyhmmer hmmsearch with a fixed Z (E-values are not used), acceleration
filters off (max mode, so short profiles are not dropped before scoring),
and a reporting floor in bits; the feature is the best domain bit score.
PSSM: pftools `pfsearch -f -k -a -r` (raw score of the best alignment for
every sequence).

PWM: lightmotif ungapped matrix of each region alignment; feature = best
window score (empty if the sequence is shorter than the matrix).
Framework: per candidate family, the mature region is predicted with that
family's full-precursor model (anatomy.method, optional lightmotif sites;
own-family holdout model for training positives), and the feature is the
smoothed log frequency of its cysteine framework among that family's
training mature peptides (own fold/group left out for training positives).

features.tsv columns: id, kind, source, role, label, tier, fold, group, then
one column per <superfamily>:<region>:<hmm|pssm|pwm> and
<superfamily>:FW:framework. Values: a score; empty = no score reported
(HMM: no domain; PWM: sequence shorter than the matrix); NA = no profile
available for this sequence (e.g. a holdout with too few sequences).
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pyhmmer
from pyhmmer import easel, plan7

from .build import safe
from .config import Config

log = logging.getLogger(__name__)

SCORE_FIELDS = ["id", "profile", "superfamily", "region", "method", "library", "score", "start", "end", "detail"]
META_FIELDS = ["id", "kind", "source", "role", "label", "tier", "fold", "group"]
NONE_LABEL = "NONE"


class ScoreError(RuntimeError):
    pass


def _s(x) -> str:
    return x if isinstance(x, str) else x.decode()


# ------------------------------------------------------------- sequences ----

def scoring_set(entries: list[dict], partition: list[dict], negatives: list[dict]) -> list[dict]:
    seqs = {e["accession"]: e["sequence"] for e in entries}
    out = []
    for p in partition:
        if p["scope"] == "conopeptide" and p["role"] in ("train", "test", "temporal_test"):
            out.append({"id": p["accession"], "kind": "positive", "source": "uniprot", "role": p["role"],
                        "label": p["superfamily"], "tier": p["tier"], "fold": p["fold"], "group": p["group"],
                        "sequence": seqs[p["accession"]]})
    for n in negatives:
        if n["role"] in ("train", "test"):
            out.append({"id": n["id"], "kind": "negative", "source": n["source"], "role": n["role"],
                        "label": NONE_LABEL, "tier": "", "fold": "", "group": n["group"], "sequence": n["sequence"]})
    return sorted(out, key=lambda r: r["id"])


def unit_of(row: dict) -> str | None:
    """Holdout unit of a training positive, as named by the build stage."""
    if row["kind"] != "positive" or row["role"] != "train":
        return None
    if row["fold"].isdigit():
        return f"fold{row['fold']}"
    if row["fold"] == "logo":
        return f"group-{safe(row['group'])}"
    return None


# -------------------------------------------------------------------- hmm ----

def hmm_scores(hmm_path: Path, seqs: list[dict], library: str, cfg: Config) -> list[dict]:
    if not seqs or not hmm_path.exists() or hmm_path.stat().st_size == 0:
        return []
    sc = cfg.score
    alpha = easel.Alphabet.amino()
    with plan7.HMMFile(hmm_path) as fh:
        hmms = list(fh)
    block = easel.TextSequenceBlock(
        [easel.TextSequence(name=r["id"].encode(), sequence=_clean(r["sequence"])) for r in seqs]).digitize(alpha)
    kw = dict(Z=1, domZ=1, T=sc.hmm_min_bits, domT=sc.hmm_min_bits, cpus=sc.threads or os.cpu_count() or 1)
    if sc.max_mode:
        kw |= dict(F1=1.0, F2=1.0, F3=1.0, bias_filter=False)
    out = []
    for hits in pyhmmer.hmmsearch(hmms, block, **kw):
        name = _s(hits.query.name)
        _, sf, region = _split_name(name)
        for h in hits:
            doms = [d for d in h.domains if d.score >= sc.hmm_min_bits]
            if not doms:
                continue
            d = max(doms, key=lambda d: (d.score, -d.env_from))
            out.append({"id": _s(h.name), "profile": name, "superfamily": sf, "region": region, "method": "hmm",
                        "library": library, "score": round(d.score, 2), "start": d.env_from, "end": d.env_to})
    return out


def _clean(seq: str) -> str:
    return re.sub(r"[^A-Z]", "X", seq.upper())


def _split_name(name: str) -> tuple[str, str, str]:
    m = re.match(r"^(CONOPEP)_(.+)_(SIG|PRO|MAT|POST)$", name)
    if not m:
        raise ScoreError(f"Unexpected profile name {name!r}")
    return m.group(1), m.group(2), m.group(3)


# ------------------------------------------------------------------- pssm ----

# xPSA header lines, with or without the "/start-end" location, e.g.
#   >q0 motif=CONOPEP_A_MAT|CONOPEP_A_MAT raw_score=911
#   >q0/12-40 level=NA raw_score=911 seq_end=-5
_PSA = re.compile(r"^>(\S+?)(?:/(-?\d+)-(-?\d+))?\s.*?raw_score=\s*(-?\d+)")


def pssm_scores(prf: Path, name: str, seqs: list[dict], library: str, pfsearch: str, options: list[str]) -> list[dict]:
    if not seqs:
        return []
    _, sf, region = _split_name(name)
    with tempfile.TemporaryDirectory(prefix="conodb-ps-") as tmp:
        fa = Path(tmp) / "seqs.fa"
        # short ids: pftools truncates long identifiers
        alias = {f"q{i}": r["id"] for i, r in enumerate(seqs)}
        fa.write_text("".join(f">q{i}\n{_clean(r['sequence'])}\n" for i, r in enumerate(seqs)))
        cmd = [pfsearch, "-f", "-k", "-a", "-r", *options, str(prf), str(fa)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise ScoreError(f"{' '.join(cmd)} failed (exit {res.returncode}): {res.stderr.strip()[:500]}")
    best: dict[str, tuple[int, int | str, int | str]] = {}
    for line in res.stdout.splitlines():
        m = _PSA.match(line)
        if not m:
            continue
        sid, sc = alias.get(m.group(1)), int(m.group(4))
        b = int(m.group(2)) if m.group(2) else ""
        e = int(m.group(3)) if m.group(3) else ""
        if sid is not None and (sid not in best or sc > best[sid][0]):
            best[sid] = (sc, b, e)
    if len(best) != len(seqs):
        sample = "\n".join(res.stdout.splitlines()[:3])
        raise ScoreError(
            f"pfsearch reported {len(best)} of {len(seqs)} sequences for {name}; output format not understood. "
            f"First lines:\n{sample}")
    return [{"id": sid, "profile": name, "superfamily": sf, "region": region, "method": "pssm", "library": library,
             "score": s, "start": b, "end": e} for sid, (s, b, e) in best.items()]


# ------------------------------------------------------------------- main ----

def pwm_scores(lib_path: Path, seqs: list[dict], library: str, pseudocount: float) -> list[dict]:
    """Best ungapped window score of each sequence for each lightmotif matrix."""
    import json

    import lightmotif

    if not seqs or not lib_path.exists():
        return []
    with open(lib_path) as fh:
        lib = json.load(fh)
    out = []
    for name, windows in sorted(lib.items()):
        if not windows:
            continue
        _, sf, region = _split_name(name)
        pssm = lightmotif.create(windows, protein=True).counts.normalize(pseudocount).log_odds()
        L = len(windows[0])
        for r in seqs:
            q = re.sub("[^ACDEFGHIKLMNPQRSTVWY]", "A", r["sequence"].upper())
            if len(q) < L:
                continue
            sc = pssm.calculate(lightmotif.stripe(q, protein=True))
            vals = [float(sc[i]) for i in range(len(q) - L + 1)]
            k = max(range(len(vals)), key=lambda i: vals[i])
            out.append({"id": r["id"], "profile": name, "superfamily": sf, "region": region, "method": "pwm",
                        "library": library, "score": round(vals[k], 3), "start": k + 1, "end": k + L})
    return out


def load_full_library(stem: Path) -> tuple[list, dict[str, list[str]]]:
    import json

    from pyhmmer import plan7
    hp, lp = Path(f"{stem}.full.hmm"), Path(f"{stem}.full.labels.json")
    if not hp.exists() or hp.stat().st_size == 0 or not lp.exists():
        return [], {}
    with plan7.HMMFile(hp) as fh:
        hmms = list(fh)
    with open(lp) as fh:
        labels = json.load(fh)
    return hmms, labels


def load_sites(path: Path, pseudocount: float):
    import json

    from .core import SiteModels
    if not path.exists():
        return None
    with open(path) as fh:
        d = json.load(fh)
    return SiteModels.from_windows(d["windows"], d["left"], d["right"], pseudocount)


def framework_scores(cfg: Config, rdir: Path, seqs: list[dict], partition: list[dict], regions: list[dict],
                     families: list[str], region_idx: dict) -> list[dict]:
    """Framework feature rows (method 'framework', region 'FW') for every sequence x family."""
    from collections import Counter

    from .core import framework_logprob, mature_framework, predict_segments

    an, sc = cfg.anatomy, cfg.score
    part = {p["accession"]: p for p in partition}
    sites = load_sites(rdir / "profiles" / "final.sites.json", an.pseudocount) if an.site_types else None
    kw = an.region_kwargs()

    # framework counts of training mature peptides, per family and per unit
    counts: dict[str, Counter] = defaultdict(Counter)
    unit_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    vocab = set()
    for r in regions:
        p = part.get(r["accession"])
        if r["region"] != "MAT" or not p or p["role"] != "train" or p["scope"] != "conopeptide":
            continue
        fw = r["framework"] or "-"
        u = unit_of({"kind": "positive", "role": "train", "fold": p["fold"], "group": p["group"]})
        counts[safe(p["superfamily"])][fw] += 1
        unit_counts[(safe(p["superfamily"]), u)][fw] += 1
        vocab.add(fw)
    V = len(vocab) + 1

    def segments(r: dict, fam: str, library: str, full_hits: dict, labels: dict) -> dict:
        name = f"CONOPEP_{fam}_FULL"
        if an.method == "transfer":
            return predict_segments(r["sequence"], full_hit=full_hits.get((r["id"], name)),
                                    labels=labels.get(name), sites=sites, **kw)
        rh = {rg: region_idx[(r["id"], fam, rg, "hmm", library)][1:]
              for rg in ("SIG", "PRO", "MAT", "POST") if (r["id"], fam, rg, "hmm", library) in region_idx
              and region_idx[(r["id"], fam, rg, "hmm", library)][0] > 0}
        return predict_segments(r["sequence"], region_hits=rh, sites=sites, **kw)

    out = []
    from .core import best_domains
    final_hmms, final_labels = load_full_library(rdir / "profiles" / "final")
    final_hits = best_domains(final_hmms, ((r["id"], r["sequence"]) for r in seqs), min_bits=sc.hmm_min_bits,
                              max_mode=sc.max_mode, cpus=sc.threads, alignments=True) if an.method == "transfer" else {}
    by_unit: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in seqs:
        u = unit_of(r)
        if u:
            by_unit[(safe(r["label"]), u)].append(r)
    unit_hits: dict[tuple[str, str], tuple[dict, dict]] = {}
    if an.method == "transfer":
        for (fam, u), rs in by_unit.items():
            hmms, labels = load_full_library(rdir / "profiles" / "holdout" / fam / u)
            hits = best_domains(hmms, ((r["id"], r["sequence"]) for r in rs), min_bits=sc.hmm_min_bits,
                                max_mode=sc.max_mode, cpus=sc.threads, alignments=True)
            unit_hits[(fam, u)] = (hits, labels)
    for r in seqs:
        u = unit_of(r)
        for fam in families:
            own = u is not None and safe(r["label"]) == fam
            if own:
                library = f"holdout:{u}"
                hits, labels = unit_hits.get((fam, u), ({}, {}))
                c = counts[fam] - unit_counts[(fam, u)]
            else:
                library, hits, labels, c = "final", final_hits, final_labels, counts[fam]
            seg = segments(r, fam, library, hits, labels)
            fw = mature_framework(r["sequence"], seg, an.short_peptide_length)
            mat = seg.get("MAT", ("", ""))
            out.append({"id": r["id"], "profile": f"CONOPEP_{fam}_FW", "superfamily": fam, "region": "FW",
                        "method": "framework", "library": library,
                        "score": round(framework_logprob(fw, c, V, sc.framework_pseudocount), 4),
                        "start": mat[0], "end": mat[1], "detail": "" if fw is None else fw})
    return out


def run_score(cfg: Config, rdir: Path, entries: list[dict], partition: list[dict], negatives: list[dict],
              profiles: list[dict], regions: list[dict] | None = None) -> tuple[list[dict], list[dict], list[str], dict]:
    sc = cfg.score
    seqs = scoring_set(entries, partition, negatives)
    built = [p for p in profiles if p["status"] == "built" and p["region"] in ("SIG", "PRO", "MAT", "POST")
             and p["set"] in ("final", "holdout")]
    use_pssm = bool(built) and all(p["pssm"] for p in built) and sc.pssm != "false"
    pfsearch = shutil.which(sc.pfsearch) if use_pssm else None
    if use_pssm and pfsearch is None:
        if sc.pssm == "true":
            raise ScoreError(f"score.pssm = true but {sc.pfsearch!r} is not on PATH")
        log.warning("PSSMs were built but %r is not on PATH: scoring HMMs only", sc.pfsearch)
        use_pssm = False
    use_pwm = sc.pwm and bool(built) and all(p.get("pwm") for p in built)

    scores: list[dict] = []
    # 1. final libraries: everyone
    scores += hmm_scores(rdir / "profiles" / "final.hmm", seqs, "final", cfg)
    # 2. holdouts: training positives against their own family's unit
    by_unit: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in seqs:
        u = unit_of(r)
        if u:
            by_unit[(r["label"], u)].append(r)
    for (sf, unit), rs in sorted(by_unit.items()):
        lib = rdir / "profiles" / "holdout" / safe(sf) / f"{unit}.hmm"
        scores += hmm_scores(lib, rs, f"holdout:{unit}", cfg)

    if use_pssm:
        jobs = [(rdir / p["pssm"], p["name"], seqs, "final") for p in built if p["set"] == "final"]
        for (sf, unit), rs in sorted(by_unit.items()):
            jobs += [(rdir / p["pssm"], p["name"], rs, f"holdout:{unit}") for p in built
                     if p["set"] == "holdout" and p["superfamily"] == sf and p["unit"] == unit]
        with ThreadPoolExecutor(max_workers=sc.threads or os.cpu_count() or 1) as pool:
            for res in pool.map(lambda j: pssm_scores(*j, pfsearch, sc.pfsearch_options), jobs):
                scores += res
    if use_pwm:
        pc = cfg.anatomy.pseudocount
        scores += pwm_scores(rdir / "profiles" / "final.pwm.json", seqs, "final", pc)
        for (sf, unit), rs in sorted(by_unit.items()):
            scores += pwm_scores(rdir / "profiles" / "holdout" / safe(sf) / f"{unit}.pwm.json", rs,
                                 f"holdout:{unit}", pc)

    # 3. feature table
    methods = ["hmm"] + (["pssm"] if use_pssm else []) + (["pwm"] if use_pwm else [])
    finals = sorted({(p["superfamily"], p["region"]) for p in built if p["set"] == "final"},
                    key=lambda k: (k[0], ["SIG", "PRO", "MAT", "POST"].index(k[1])))
    families = sorted({safe(sf) for sf, _ in finals})
    holdout_have = {(p["superfamily"], p["unit"], p["region"]) for p in built if p["set"] == "holdout"}
    idx: dict[tuple[str, str, str, str, str], tuple] = {}
    for s_ in scores:
        idx[(s_["id"], safe(s_["superfamily"]), s_["region"], s_["method"], s_["library"])] = (
            s_["score"], s_["start"], s_["end"])
    use_fw = sc.framework and regions is not None and families
    if use_fw:
        fw_rows = framework_scores(cfg, rdir, seqs, partition, regions, families, idx)
        scores += fw_rows
        for s_ in fw_rows:
            idx[(s_["id"], s_["superfamily"], "FW", "framework", s_["library"])] = (s_["score"], s_["start"], s_["end"])
    cols = [f"{safe(sf)}:{rg}:{m}" for sf, rg in finals for m in methods]
    if use_fw:
        cols += [f"{fam}:FW:framework" for fam in families]
    features = []
    n_na = 0
    for r in seqs:
        row = {k: r[k] for k in META_FIELDS}
        u = unit_of(r)
        for sf, rg in finals:
            own = u is not None and sf == r["label"]
            for m in methods:
                col = f"{safe(sf)}:{rg}:{m}"
                if own:
                    if (sf, u, rg) not in holdout_have:
                        row[col] = "NA"
                        n_na += 1
                        continue
                    lib = f"holdout:{u}"
                else:
                    lib = "final"
                v = idx.get((r["id"], safe(sf), rg, m, lib))
                row[col] = "" if v is None else v[0]
        if use_fw:
            for fam in families:
                lib = f"holdout:{u}" if (u is not None and safe(r["label"]) == fam) else "final"
                v = idx.get((r["id"], fam, "FW", "framework", lib))
                row[f"{fam}:FW:framework"] = "" if v is None else v[0]
        features.append(row)
    scores.sort(key=lambda s_: (s_["id"], s_["profile"], s_["method"], s_["library"]))
    stats = {"sequences": len(seqs), "positives": sum(r["kind"] == "positive" for r in seqs),
             "negatives": sum(r["kind"] == "negative" for r in seqs), "features": len(cols),
             "scores": len(scores), "pssm": use_pssm, "pwm": use_pwm, "framework": bool(use_fw),
             "no_profile_cells": n_na}
    return scores, features, META_FIELDS + cols, stats
