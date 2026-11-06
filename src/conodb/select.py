"""Select: choose the defaults by cross-validation on training rows only.

The test set is evaluated once, with defaults fixed beforehand; this stage is
where they are chosen. Everything here uses training rows (role 'train') and
features that are already out-of-fold with respect to the profiles.

Classifier variants (`variants.classifier_variants`)
  Grouped K-fold cross-validation: tier-1 positives keep their partition
  fold; tier-2 entry groups, real negatives and their groups are assigned to
  folds by a stable hash of the group; a shuffled negative goes with the fold
  of the precursor it was made from. For each variant, fold and tau in
  evaluate.select_taus, the penalised model is refitted on the other folds
  (its own standardisation and fills) and predicts the held-out fold.
  Primary metric: log loss, weighted as in the fit (1 / entry-group size,
  shuffled and real negatives balanced).

Region variants (`variants.region_variants`)
  Every training conopeptide is cut with its own family's holdout models
  (the full-precursor and region HMMs built without its fold or group) and
  site matrices rebuilt without the windows of that same unit. Primary
  metric: boundary_exact (signal end, mature start and mature end predicted
  exactly, pooled over the EXACT UniProt boundaries).

Intervals: percentile bootstrap over entry groups (evaluate.bootstrap
resamples), paired between variants.

Rule (written down before looking at the test set): a variant is eligible
when its paired difference to the best variant has a 95% interval that
includes zero (not clearly worse), overall and on sequences without a signal
peptide (mature peptides and fragments: a minority that an overall metric
could trade away; this guard applies when that subset has >= 10 entry
groups). Among eligible variants the simplest (`complexity`) wins, ties
broken by the overall point estimate. Classifier variants are chosen at tau = fit.tau_fixed; the
other taus show whether the choice depends on shrinkage.
"""

from __future__ import annotations

import hashlib
import logging
import multiprocessing as mp
import os
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np

from .anatomy import aggregate, score_sequence, truth_of
from .build import safe
from .config import Config
from .core import SiteModels, best_domains, predict_segments
from .fit import NothingToFit, fit_map, make_design, predict_proba
from .io import read_tsv
from .score import load_full_library, unit_of
from .variants import classifier_variants, region_kwargs, region_variants

log = logging.getLogger(__name__)

MIN_GUARD_GROUPS = 10        # entry groups needed in a subset for it to guard the choice
CLS_FIELDS = ["variant", "note", "methods", "regions", "hmm_floor", "tau", "complexity", "n", "log_loss",
              "log_loss_lo", "log_loss_hi", "diff_vs_best_lo", "diff_vs_best_hi", "no_signal_log_loss",
              "no_signal_diff_lo", "no_signal_diff_hi", "not_worse", "selected",
              "accuracy", "positives_correct", "positives_correct_lo", "positives_correct_hi", "family_if_known",
              "negatives_called_family", "macro_recall", "complete_correct", "no_signal_correct",
              "no_signal_family_if_known"]
ERR_FIELDS = ["variant", "boundary", "direction", "residues", "n"]
REG_FIELDS = ["variant", "note", "method", "extension", "overhang", "sites", "cterm", "complexity", "n", "boundary_exact",
              "boundary_exact_lo", "boundary_exact_hi", "diff_vs_best_lo", "diff_vs_best_hi",
              "no_signal_boundary_exact", "no_signal_diff_lo", "no_signal_diff_hi", "not_worse", "selected",
              "boundary_n", "signal_end_exact", "mature_start_exact", "mature_end_exact", "mature_exact",
              "mature_exact_n", "residue_accuracy", "complete_mature_exact", "no_signal_mature_exact"]


def _hash_fold(key: str, k: int) -> int:
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16) % k + 1


def cv_folds(rows: list[dict], k: int) -> dict[str, int]:
    """Fold (1..k) of every training row; groups never straddle folds."""
    group_fold = {}
    for r in rows:
        if r["kind"] == "positive" and r["fold"].isdigit():
            group_fold[r["group"]] = int(r["fold"])
    out = {}
    for r in rows:
        if r["kind"] == "positive":
            out[r["id"]] = int(r["fold"]) if r["fold"].isdigit() else group_fold.setdefault(
                r["group"], _hash_fold("pos:" + r["group"], k))
        elif r["group"].startswith("SHUF-") and r["group"][5:] in group_fold:
            out[r["id"]] = group_fold[r["group"][5:]]
        else:
            out[r["id"]] = _hash_fold(f"{r['kind']}:{r['group']}", k)
    return out


def bootstrap_counts(n_groups: int, B: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.multinomial(n_groups, np.full(n_groups, 1 / n_groups), size=B).astype(float) if B else \
        np.zeros((0, n_groups))


def boot_ratio(num: np.ndarray, den: np.ndarray, gidx: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Bootstrap distribution (B x C) of sum(num[:, c]) / sum(den) with groups resampled."""
    G = counts.shape[1]
    ng = np.zeros((G, num.shape[1]))
    np.add.at(ng, gidx, num)
    dg = np.zeros(G)
    np.add.at(dg, gidx, den)
    d = counts @ dg
    with np.errstate(invalid="ignore", divide="ignore"):
        return (counts @ ng) / d[:, None]


def _ci(x: np.ndarray) -> tuple[float, float]:
    x = x[~np.isnan(x)]
    if not x.size:
        return (np.nan, np.nan)
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


def vs_best(point: np.ndarray, boot: np.ndarray, higher_better: bool) -> tuple[int, np.ndarray, np.ndarray]:
    """Best index, paired-difference interval (2 x C) to it, and 'not clearly worse'."""
    sign = 1 if higher_better else -1
    best = int(np.nanargmax(sign * point))
    C = len(point)
    lo, hi = np.full(C, np.nan), np.full(C, np.nan)
    if boot.size:
        for c in range(C):
            lo[c], hi[c] = _ci(boot[:, c] - boot[:, best])
        nw = (hi >= 0) if higher_better else (lo <= 0)
    else:
        nw = np.arange(C) == best
    nw = nw & ~np.isnan(point)
    nw[best] = True
    return best, np.vstack([lo, hi]), nw


def choose(point: np.ndarray, boot: np.ndarray, complexity: list[int], higher_better: bool,
           guards: tuple = ()) -> tuple[int, int, np.ndarray, np.ndarray]:
    """Best index (primary metric), selected index, primary difference interval, and
    eligibility: not clearly worse than the best on the primary metric and on
    every guard (point, boot) metric (e.g. the same metric on mature peptides
    alone, so a minority subset cannot be traded away)."""
    sign = 1 if higher_better else -1
    best, diff, eligible = vs_best(point, boot, higher_better)
    for gp, gb in guards:
        eligible = eligible & vs_best(gp, gb, higher_better)[2]
    if not eligible.any():
        eligible[best] = True
    cand = [c for c in range(len(point)) if eligible[c]]
    sel = min(cand, key=lambda c: (complexity[c], -sign * point[c]))
    return best, sel, diff, eligible


# ------------------------------------------------------------- classifier ----

_G: dict = {}
_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def _init_worker(g: dict) -> None:
    _G.update(g)


def run_parallel(func, tasks: list, shared: dict, n_proc: int) -> list:
    """Map func over tasks in fresh (spawned) processes, one BLAS thread each;
    `shared` is sent once per process. Serial when n_proc <= 1."""
    if n_proc <= 1:
        _G.update(shared)
        try:
            return [func(t) for t in tasks]
        finally:
            _G.clear()
    saved = {k: os.environ.get(k) for k in _THREAD_VARS}
    os.environ.update({k: "1" for k in _THREAD_VARS})
    try:
        with mp.get_context("spawn").Pool(n_proc, initializer=_init_worker, initargs=(shared,)) as pool:
            return pool.map(func, tasks, chunksize=1)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _cls_task(args):
    vi, fold = args
    g = _G
    v = g["variants"][vi]
    cfg = g["cfg"]
    cv = replace(cfg, fit=replace(cfg.fit, methods=list(v["methods"]), regions=list(v["regions"]),
                                  hmm_floor=v["hmm_floor"]))
    rows = [r | {"role": "train" if g["fold"][r["id"]] != fold else "test"} for r in g["rows"]]
    try:
        d = make_design(cv, rows, g["fields"])
    except NothingToFit:
        return vi, fold, None
    test = np.array([ro == "test" for ro in d.roles])
    out = {}
    for tau in g["taus"]:
        m = fit_map(cv, d, tau)
        out[tau] = predict_proba(d, m["alpha"], m["B"])[test]
    return vi, fold, ([i for i, t in zip(d.ids, test, strict=True) if t], out)


def _pool(cores: int) -> int:
    return cores if cores > 0 else (os.cpu_count() or 1)


def select_classifier(cfg: Config, rdir: Path, types_of: dict[str, str]) -> dict:
    ec, fc = cfg.evaluate, cfg.fit
    rows = [r for r in read_tsv(rdir / "features.tsv") if r["role"] == "train"]
    with open(rdir / "features.tsv", encoding="utf-8") as fh:
        fields = fh.readline().rstrip("\n").split("\t")
    have = {c.split(":")[2] for c in fields if c.count(":") == 2}
    variants = classifier_variants(cfg, have)
    taus = sorted(set(ec.select_taus) | {fc.tau_fixed})
    k = cfg.partition.folds
    folds = cv_folds(rows, k)
    # weights and classes from one design with every training row fitted
    d0 = make_design(cfg, rows, fields)
    w = dict(zip(d0.ids, d0.w, strict=True))
    K = len(d0.families)
    yid = dict(zip(d0.ids, d0.y, strict=True))
    ids = d0.ids
    pos = {i: n for n, i in enumerate(ids)}

    shared = {"cfg": cfg, "rows": rows, "fields": fields, "variants": variants, "fold": folds, "taus": taus}
    tasks = [(vi, f) for vi in range(len(variants)) for f in range(1, k + 1)]
    n_proc = min(_pool(fc.cores), len(tasks))
    log.info("  classifier: %d variants x %d folds x %d taus (%d processes)", len(variants), k, len(taus), n_proc)
    probs = {(vi, t): np.full((len(ids), K + 1), np.nan) for vi in range(len(variants)) for t in taus}
    results = run_parallel(_cls_task, tasks, shared, n_proc)
    usable = set(range(len(variants)))
    for vi, _, res in results:
        if res is None:
            usable.discard(vi)
            continue
        tid, out = res
        rows_i = [pos[i] for i in tid]
        for t, pr in out.items():
            probs[(vi, t)][rows_i] = pr
    variants_u = [v for vi, v in enumerate(variants) if vi in usable]
    vidx = [vi for vi in range(len(variants)) if vi in usable]

    y = np.array([yid[i] for i in ids])
    wv = np.array([w[i] for i in ids])
    rowmap = {r["id"]: r for r in rows}
    groups = [f"{rowmap[i]['kind']}:{rowmap[i]['group']}" for i in ids]
    gid = {g: n for n, g in enumerate(dict.fromkeys(groups))}
    gidx = np.array([gid[g] for g in groups])
    counts = bootstrap_counts(len(gid), ec.bootstrap, ec.seed)
    typ = np.array([types_of.get(i, "") for i in ids])
    posm = y < K

    table = []
    sel_info = {}
    for t in taus:
        C = len(vidx)
        loss = np.zeros((len(ids), C))
        correct = np.zeros((len(ids), C))
        for c, vi in enumerate(vidx):
            p = np.clip(probs[(vi, t)], 1e-12, 1)
            loss[:, c] = -np.log(p[np.arange(len(ids)), y])
            correct[:, c] = (p.argmax(1) == y)
        point = (wv[:, None] * loss).sum(0) / wv.sum()
        boot = boot_ratio(wv[:, None] * loss, wv, gidx, counts)
        pc_boot = boot_ratio(correct * posm[:, None], posm.astype(float), gidx, counts)
        nm_pos = (posm & (typ == "no signal peptide")).astype(float)
        guards, g_point, g_diff = (), np.full(C, np.nan), np.full((2, C), np.nan)
        if len(set(gidx[nm_pos > 0])) >= MIN_GUARD_GROUPS:
            g_point = (wv[:, None] * loss * nm_pos[:, None]).sum(0) / (wv * nm_pos).sum()
            g_boot = boot_ratio(wv[:, None] * loss * nm_pos[:, None], wv * nm_pos, gidx, counts)
            guards = ((g_point, g_boot),)
            g_diff = vs_best(g_point, g_boot, False)[1]
        best, sel, diff, not_worse = choose(point, boot, [v["complexity"] for v in variants_u], False, guards)
        sel_info[t] = {"best": variants_u[best]["name"], "selected": variants_u[sel]["name"],
                       "order": [variants_u[c]["name"] for c in np.argsort(point)]}
        for c, (vi, v) in enumerate(zip(vidx, variants_u, strict=True)):
            p = probs[(vi, t)]
            pred = p.argmax(1)
            fam = p[:, :K].argmax(1)
            rec = [float(np.mean(pred[y == cl] == cl)) for cl in sorted(set(y[posm].tolist()))]
            lo, hi = _ci(boot[:, c]) if boot.size else (np.nan, np.nan)
            plo, phi = _ci(pc_boot[:, c]) if pc_boot.size else (np.nan, np.nan)
            cm, nm = posm & (typ == "complete precursor"), posm & (typ == "no signal peptide")
            table.append({
                "variant": v["name"], "note": v["note"], "methods": "+".join(v["methods"]),
                "regions": "/".join(v["regions"]), "hmm_floor": v["hmm_floor"], "tau": t,
                "complexity": v["complexity"], "n": len(ids), "log_loss": _r(point[c]), "log_loss_lo": _r(lo),
                "log_loss_hi": _r(hi), "diff_vs_best_lo": _r(diff[0, c], 5), "diff_vs_best_hi": _r(diff[1, c], 5),
                "not_worse": bool(not_worse[c]), "selected": c == sel,
                "no_signal_log_loss": _r(g_point[c]), "no_signal_diff_lo": _r(g_diff[0, c], 5),
                "no_signal_diff_hi": _r(g_diff[1, c], 5),
                "accuracy": _r(np.mean(pred == y)), "positives_correct": _r(np.mean(pred[posm] == y[posm])),
                "positives_correct_lo": _r(plo), "positives_correct_hi": _r(phi),
                "family_if_known": _r(np.mean(fam[posm] == y[posm])),
                "negatives_called_family": _r(np.mean(pred[~posm] != K)) if (~posm).any() else "",
                "macro_recall": _r(np.mean(rec)) if rec else "",
                "complete_correct": _r(np.mean(pred[cm] == y[cm])) if cm.any() else "",
                "no_signal_correct": _r(np.mean(pred[nm] == y[nm])) if nm.any() else "",
                "no_signal_family_if_known": _r(np.mean(fam[nm] == y[nm])) if nm.any() else ""})
    # tau sensitivity: rank agreement with the selection tau, and whether the choice stays "not worse"
    ref = fc.tau_fixed
    sens = []
    ll = {(r["variant"], r["tau"]): r["log_loss"] for r in table}
    nw = {(r["variant"], r["tau"]): r["not_worse"] for r in table}
    chosen = sel_info[ref]["selected"]
    names = [v["name"] for v in variants_u]
    for t in taus:
        a = np.array([ll[(n, ref)] for n in names], dtype=float)
        b = np.array([ll[(n, t)] for n in names], dtype=float)
        sens.append({"tau": t, "best": sel_info[t]["best"], "selected_at_this_tau": sel_info[t]["selected"],
                     "choice_not_worse_here": nw[(chosen, t)],
                     "spearman_vs_selection_tau": _r(_spearman(a, b))})
    sv = next(v for v in variants_u if v["name"] == chosen)
    return {"table": table, "selected": sv, "tau": ref, "best": sel_info[ref]["best"], "sensitivity": sens,
            "n": len(ids), "folds": k, "groups": len(gid)}


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return float("nan")
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def _r(x, nd=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return round(float(x), nd)


# ----------------------------------------------------------------- regions ---

def select_regions(cfg: Config, rdir: Path) -> dict:
    import json

    an, sc, ec = cfg.anatomy, cfg.score, cfg.evaluate
    entries = {e["accession"]: e for e in read_tsv(rdir / "entries.tsv")}
    part = read_tsv(rdir / "partition.tsv")
    truth, exact = truth_of(read_tsv(rdir / "regions.tsv"))
    variants = region_variants(cfg)
    site_path = rdir / "profiles" / "final.sites.json"
    sites_doc = json.loads(site_path.read_text()) if site_path.exists() else None
    if sites_doc is not None and "accessions" not in sites_doc:
        raise RuntimeError("profiles/final.sites.json has no window accessions (built by an older conodb): "
                           "rerun with --from build")

    units: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for p in part:
        if p["scope"] != "conopeptide" or p["role"] != "train" or p["tier"] not in ("1", "2"):
            continue
        u = unit_of({"kind": "positive", "role": "train", "fold": p["fold"], "group": p["group"]})
        if u:
            units[(safe(p["superfamily"]), u)].append(p)
    from pyhmmer import plan7

    items: dict[str, list[dict]] = defaultdict(list)
    for (fam, u), ps in sorted(units.items()):
        stem = rdir / "profiles" / "holdout" / fam / u
        full_hmms, labels = load_full_library(stem)
        hmm_path = Path(f"{stem}.hmm")
        reg_hmms = []
        if hmm_path.exists() and hmm_path.stat().st_size:
            with plan7.HMMFile(hmm_path) as fh:
                reg_hmms = list(fh)
        if not full_hmms and not reg_hmms:
            continue
        seqs = [(p["accession"], entries[p["accession"]]["sequence"]) for p in ps]
        hf = best_domains(full_hmms, seqs, min_bits=sc.hmm_min_bits, max_mode=sc.max_mode, cpus=sc.threads,
                          alignments=True)
        hr = best_domains(reg_hmms, seqs, min_bits=sc.hmm_min_bits, max_mode=sc.max_mode, cpus=sc.threads)
        sites = None
        if sites_doc:
            own = {p["accession"] for p in ps}
            win = {t: [w for w, a in zip(ws, sites_doc["accessions"].get(t, []), strict=True) if a not in own]
                   for t, ws in sites_doc["windows"].items()}
            sites = SiteModels.from_windows(win, sites_doc["left"], sites_doc["right"], an.pseudocount)
        name = f"CONOPEP_{fam}_FULL"
        for p in ps:
            acc = p["accession"]
            seq = entries[acc]["sequence"]
            env = {rg: (h.start, h.end) for rg in ("SIG", "PRO", "MAT", "POST")
                   if (h := hr.get((acc, f"CONOPEP_{fam}_{rg}"))) is not None and h.score > 0}
            typ = "complete precursor" if entries[acc]["complete_precursor"] == "true" else (
                "signal, incomplete" if "SIG" in truth[acc] else "no signal peptide")
            for v in variants:
                seg = predict_segments(seq, full_hit=hf.get((acc, name)), labels=labels.get(name), region_hits=env,
                                       sites=sites if v["sites"] else None, **region_kwargs(cfg, v))
                it = score_sequence(acc, seq, seg, truth[acc], exact)
                it.update({"group": p["group"], "type": typ, "seq": seq})
                items[v["name"]].append(it)

    names = [v["name"] for v in variants if items.get(v["name"])]
    if not names:
        return {"table": [], "selected": None, "n": 0}
    base = items[names[0]]
    gid = {g: n for n, g in enumerate(dict.fromkeys(it["group"] for it in base))}
    gidx = np.array([gid[it["group"]] for it in base])
    counts = bootstrap_counts(len(gid), ec.bootstrap, ec.seed + 1)
    hits = np.zeros((len(base), len(names)))
    for c, n in enumerate(names):
        hits[:, c] = [sum(it["boundaries"].values()) for it in items[n]]
    den = np.array([len(it["boundaries"]) for it in base], dtype=float)
    point = hits.sum(0) / den.sum() if den.sum() else np.full(len(names), np.nan)
    boot = boot_ratio(hits, den, gidx, counts)
    vmap = {v["name"]: v for v in variants}
    nos = np.array([it["type"] == "no signal peptide" for it in base], dtype=float)
    guards, g_point, g_diff = (), np.full(len(names), np.nan), np.full((2, len(names)), np.nan)
    if len(set(gidx[nos > 0])) >= MIN_GUARD_GROUPS and (den * nos).sum():
        g_point = (hits * nos[:, None]).sum(0) / (den * nos).sum()
        g_boot = boot_ratio(hits * nos[:, None], den * nos, gidx, counts)
        guards = ((g_point, g_boot),)
        g_diff = vs_best(g_point, g_boot, True)[1]
    best, sel, diff, not_worse = choose(point, boot, [vmap[n]["complexity"] for n in names], True, guards)
    table = []
    for c, n in enumerate(names):
        v = vmap[n]
        agg = aggregate(items[n])
        lo, hi = _ci(boot[:, c]) if boot.size else (np.nan, np.nan)
        comp = aggregate([it for it in items[n] if it["type"] == "complete precursor"])
        nos = aggregate([it for it in items[n] if it["type"] == "no signal peptide"])
        table.append({"variant": n, "note": v["note"], "method": v["method"], "extension": v["extension"],
                      "overhang": v["overhang"], "sites": "+".join(v["sites"]), "cterm": v["cterm"],
                      "complexity": v["complexity"],
                      **agg, "boundary_exact": _r(point[c]), "boundary_exact_lo": _r(lo), "boundary_exact_hi": _r(hi),
                      "diff_vs_best_lo": _r(diff[0, c], 5), "diff_vs_best_hi": _r(diff[1, c], 5),
                      "not_worse": bool(not_worse[c]), "selected": c == sel,
                      "no_signal_boundary_exact": _r(g_point[c]), "no_signal_diff_lo": _r(g_diff[0, c], 5),
                      "no_signal_diff_hi": _r(g_diff[1, c], 5),
                      "complete_mature_exact": comp["mature_exact"], "no_signal_mature_exact": nos["mature_exact"]})
    from collections import Counter

    from .anatomy import error_residues
    errors = []
    for n in dict.fromkeys([names[sel], names[best]] + [v["name"] for v in variants if "current settings" in v["note"]]):
        if n not in items:
            continue
        cnt = Counter()
        for it in items[n]:
            for kind, pc, tc in it["errors"]:
                d, res = error_residues(it["seq"], kind, pc, tc)
                cnt[(kind, d, res if len(res) <= 6 else f"{len(res)} residues")] += 1
        errors += [{"variant": n, "boundary": k, "direction": d, "residues": r, "n": c}
                   for (k, d, r), c in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))]
    return {"table": table, "selected": vmap[names[sel]], "best": names[best], "n": len(base), "groups": len(gid),
            "errors": errors}


# --------------------------------------------------------------------- run ---

def run_select(cfg: Config, rdir: Path) -> dict:
    from .evaluate import sequence_types
    from .fit import load_design

    d = load_design(cfg, rdir)
    types = dict(zip(d.ids, sequence_types(read_tsv(rdir / "entries.tsv"), read_tsv(rdir / "regions.tsv"), d),
                     strict=True))
    cls = select_classifier(cfg, rdir, types)
    reg = select_regions(cfg, rdir) if cfg.build.full_models else {"table": [], "selected": None, "n": 0}
    rec = {"fit": {"methods": list(cls["selected"]["methods"]), "regions": list(cls["selected"]["regions"]),
                   "hmm_floor": cls["selected"]["hmm_floor"]},
           "classifier_variant": cls["selected"]["name"]}
    if reg["selected"]:
        v = reg["selected"]
        rec["anatomy"] = {"method": v["method"], "sites": list(v["sites"]), "cterm_trim": v["cterm"]}
        if v["method"] == "transfer":
            rec["anatomy"] |= {"extension": v["extension"], "terminal_overhang": v["overhang"]}
        rec["region_variant"] = v["name"]
    cur = {"fit": {"methods": list(cfg.fit.methods), "regions": list(cfg.fit.regions), "hmm_floor": cfg.fit.hmm_floor},
           "anatomy": {"method": cfg.anatomy.method, "sites": list(cfg.anatomy.site_types),
                       "cterm_trim": cfg.anatomy.cterm_trim,
                       "extension": cfg.anatomy.extension, "terminal_overhang": cfg.anatomy.terminal_overhang}}
    changes = [f"{s}.{k} = {_toml(v)}" for s, kv in rec.items() if isinstance(kv, dict) for k, v in kv.items()
               if cur.get(s, {}).get(k) != v]
    return {"classifier": cls, "regions": reg, "recommended": rec, "changes": changes}


def _toml(v) -> str:
    if isinstance(v, str):
        return f'"{v}"'
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list | tuple):
        return "[" + ", ".join(_toml(x) for x in v) + "]"
    return f"{v:g}" if isinstance(v, float) else str(v)


def render_selection(cfg: Config, release: str, out: dict) -> str:
    from .audit import _table

    cls, reg, ec = out["classifier"], out["regions"], cfg.evaluate
    lines = [f"# Selection of the defaults: UniProt {release}", "",
             "Cross-validation on training rows only; the test set is not used here. Rule, fixed in advance: "
             "a variant is eligible when its paired difference to the best has a 95% bootstrap interval "
             "including zero, overall and on sequences without a signal peptide (when that subset has at least "
             f"{MIN_GUARD_GROUPS} entry groups); the simplest eligible variant wins, ties broken by the overall "
             "point estimate. Intervals: "
             f"{ec.bootstrap} resamples of entry groups.", ""]
    if out["changes"]:
        lines += ["## Recommended settings (differ from the current config)", "", "```toml"]
        sec = None
        for s, kv in out["recommended"].items():
            if not isinstance(kv, dict):
                continue
            items = [(k, v) for k, v in kv.items() if f"{s}.{k} = {_toml(v)}" in out["changes"]]
            if items:
                if sec != s:
                    lines.append(f"[{s}]")
                    sec = s
                lines += [f"{k} = {_toml(v)}" for k, v in items]
        lines += ["```", "", "Set these, rerun, and the test-set numbers in `evaluation.md` / `anatomy.md` "
                  "become the ones to report.", ""]
    else:
        lines += ["## Recommended settings", "", "The current config already matches the recommendation.", ""]
    lines += [f"## Classifier ({cls['n']} training rows, {cls['groups']} entry groups, {cls['folds']} folds, "
              f"tau = {cls['tau']:g})", "",
              f"Best by log loss: **{cls['best']}**; selected: **{cls['selected']['name']}**.", ""]
    rows = [r for r in cls["table"] if r["tau"] == cls["tau"]]
    rows.sort(key=lambda r: r["log_loss"])
    lines.append(_table([{**r, "log loss [95% CI]": f"{r['log_loss']} [{r['log_loss_lo']}, {r['log_loss_hi']}]",
                          "diff vs best": f"[{r['diff_vs_best_lo']}, {r['diff_vs_best_hi']}]",
                          "no signal: log loss, diff vs best": "" if r["no_signal_log_loss"] == "" else
                          f"{r['no_signal_log_loss']} [{r['no_signal_diff_lo']}, {r['no_signal_diff_hi']}]",
                          "pick": "**selected**" if r["selected"] else ("eligible" if r["not_worse"] else "")}
                         for r in rows],
                        ["variant", "note", "complexity", "log loss [95% CI]", "diff vs best",
                         "no signal: log loss, diff vs best", "pick",
                         "positives_correct", "family_if_known", "negatives_called_family", "complete_correct",
                         "no_signal_correct", "no_signal_family_if_known"]))
    lines += ["### Shrinkage sensitivity", "",
              "Same grid refitted at other tau values: if the selected variant stays eligible and the ranking "
              "agrees, the choice does not depend on shrinkage (so the Bayesian posterior is not needed to make it).",
              "", _table(cls["sensitivity"], ["tau", "best", "selected_at_this_tau", "choice_not_worse_here",
                                              "spearman_vs_selection_tau"])]
    if reg["table"]:
        lines += [f"## Region prediction ({reg['n']} training conopeptides, {reg['groups']} entry groups; "
                  "own-family holdout models and site matrices)", "",
                  f"Best by boundary_exact: **{reg['best']}**; selected: **{reg['selected']['name']}**.", ""]
        rr = sorted(reg["table"], key=lambda r: -(r["boundary_exact"] or 0))
        lines.append(_table([{**r, "boundary exact [95% CI]":
                              f"{r['boundary_exact']} [{r['boundary_exact_lo']}, {r['boundary_exact_hi']}]",
                              "diff vs best": f"[{r['diff_vs_best_lo']}, {r['diff_vs_best_hi']}]",
                              "no signal: boundary exact, diff vs best": "" if r["no_signal_boundary_exact"] == "" else
                              f"{r['no_signal_boundary_exact']} [{r['no_signal_diff_lo']}, {r['no_signal_diff_hi']}]",
                              "pick": "**selected**" if r["selected"] else ("eligible" if r["not_worse"] else "")}
                             for r in rr],
                            ["variant", "note", "complexity", "boundary exact [95% CI]", "diff vs best",
                             "no signal: boundary exact, diff vs best", "pick",
                             "signal_end_exact", "mature_start_exact", "mature_end_exact", "mature_exact",
                             "residue_accuracy", "complete_mature_exact", "no_signal_mature_exact"]))
    err = [r for r in reg.get("errors", []) if reg.get("selected") and r["variant"] == reg["selected"]["name"]]
    if err:
        lines += [f"### Boundary errors of the selected variant (training cross-validation; top 15 of {len(err)})", "",
                  "`residues` lie between the predicted and the true boundary; `extra` = kept in the predicted "
                  "mature peptide although outside the true one.", "",
                  _table(err[:15], ["boundary", "direction", "residues", "n"])]
    lines += ["All rows (every tau, every metric): `selection_classification.tsv`, `selection_regions.tsv`, "
              "`selection_region_errors.tsv`; "
              "the recommendation: `selection_recommended.json`."]
    return "\n".join(lines) + "\n"
