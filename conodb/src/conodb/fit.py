"""Fit: baseline cutoffs, penalised and Bayesian multinomial classifiers.

Input: features.tsv (out-of-fold scores, from the score stage).

Design: for sequence i and candidate family k, X[i, k, p] is the score on
family k's profile p (p = region x method), standardised per p over the
training rows. Classes: the K families with final profiles, plus NONE.

Three classifiers, all trained on 'train' rows and evaluated on 'test' and
'temporal_test' rows:

  baseline   per-profile HMM cutoffs = the (1 - fit.baseline_fpr) quantile of
             training negatives' scores; family score = sum of excesses over
             the cutoffs; highest wins, NONE if no profile passes.
  map        the conditional-logit model below, penalised maximum likelihood
             with fixed shrinkage (fit.tau_fixed): numpy + scipy, seconds.
  bayes      the same model with tau, sigma_alpha estimated (hierarchical),
             full posterior with Stan (cmdstanpy); predictions average the
             posterior. Optional: fit.bayes = "auto" | "true" | "false".

Model (see stan/conodictor_mnl.stan):
  eta[i, k] = alpha[k] + sum_p X[i, k, p] B[k, p] + sum_p M[i, k, p] gamma[p]
  B[k, p] = beta[p] + delta[k, p],  delta[k, p] ~ N(0, tau[p]),  eta[i, NONE] = 0
Case weights: 1 / size of the entry group (training rows).
"""

from __future__ import annotations

import logging
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax

from .build import safe
from .config import Config
from .io import read_tsv, write_json

log = logging.getLogger(__name__)

# Bump when the fitting code changes results for identical inputs (the Stan
# file is hashed separately), so existing releases are refitted.
FIT_FORMAT = 6

REGIONS = ("SIG", "PRO", "MAT", "POST")
NONE = "NONE"
STAN_FILE = Path(__file__).parent / "stan" / "conodictor_mnl.stan"


class FitError(RuntimeError):
    pass


class NothingToFit(FitError):
    """No profiles or no training rows (e.g. everything is held out)."""


# ----------------------------------------------------------------- design ----

@dataclass
class Design:
    ids: list[str]
    roles: list[str]
    kinds: list[str]
    sources: list[str]
    labels: list[str]
    families: list[str]          # safe names, class order 0..K-1; NONE = K
    feats: list[tuple[str, str]]  # (region, method), p order
    raw: np.ndarray               # N x K x P raw scores after floor/fill (nan = no profile / NA)
    X: np.ndarray                 # standardised, 0 where unavailable
    M: np.ndarray                 # 1 = out-of-fold profile missing (NA)
    y: np.ndarray                 # 0..K
    w: np.ndarray                 # case weights (fitting rows), else 0
    center: np.ndarray            # P
    scale: np.ndarray             # P
    available: np.ndarray | None = None   # K x P: family k has feature p
    fill: np.ndarray | None = None        # P: value for "no score" (non-HMM features)
    hmm_floor: float = 0.0

    def feature_names(self) -> list[str]:
        return [f"{r}:{m}" for r, m in self.feats]


FEATURE_ORDER = {"SIG": 0, "PRO": 1, "MAT": 2, "POST": 3, "FW": 4}


def make_design(cfg: Config, rows: list[dict], fields: list[str], fit_roles: tuple[str, ...] = ("train",)) -> Design:
    """Feature tensor from features.tsv.

    Same transformation as core.apply_model at prediction time: HMM scores
    below fit.hmm_floor (or not reported) are set to it; other missing
    scores get the lowest value seen on the fitting rows; each feature is
    standardised on the fitting rows; a feature a family does not have is 0.
    """
    fc = cfg.fit
    methods = fc.methods
    cols = [c for c in fields if c.count(":") == 2]
    fams = sorted({c.split(":")[0] for c in cols})
    feats = sorted({(c.split(":")[1], c.split(":")[2]) for c in cols
                    if c.split(":")[2] in methods and (c.split(":")[1] in fc.regions or c.split(":")[1] == "FW")},
                   key=lambda f: (FEATURE_ORDER[f[0]], methods.index(f[1])))
    if not fams or not feats:
        raise NothingToFit("features.tsv has no usable feature columns (no profiles were built)")
    rows = [r for r in rows if r["role"] in ("train", "test", "temporal_test")]
    fam_idx = {f: k for k, f in enumerate(fams)}
    keep, skipped = [], Counter()
    for r in rows:
        if r["label"] != NONE and safe(r["label"]) not in fam_idx:
            skipped[r["label"]] += 1
            continue
        keep.append(r)
    if skipped:
        log.warning("Positives of families without profiles left out: %s", dict(skipped))
    N, K, P = len(keep), len(fams), len(feats)
    raw = np.full((N, K, P), np.nan)
    below = np.zeros((N, K, P), dtype=bool)
    M = np.zeros((N, K, P))
    available = np.zeros((K, P), dtype=bool)
    present = set(fields)
    for k, f in enumerate(fams):
        for p, (rg, m) in enumerate(feats):
            available[k, p] = f"{f}:{rg}:{m}" in present
    for i, r in enumerate(keep):
        for k, f in enumerate(fams):
            for p, (rg, m) in enumerate(feats):
                if not available[k, p]:
                    continue
                v = r.get(f"{f}:{rg}:{m}")
                if v == "NA":            # holdout profile could not be built
                    M[i, k, p] = 1
                elif v == "" or v is None:
                    below[i, k, p] = True
                else:
                    raw[i, k, p] = float(v)
    fitrows = np.array([r["role"] in fit_roles for r in keep])
    fill = np.full(P, np.nan)
    for p, (_, m) in enumerate(feats):
        if m == "hmm":
            raw[:, :, p] = np.where(below[:, :, p], fc.hmm_floor, np.maximum(raw[:, :, p], fc.hmm_floor))
        else:
            vals = raw[fitrows, :, p]
            vals = vals[~np.isnan(vals)]
            fill[p] = float(vals.min()) if vals.size else 0.0
            raw[:, :, p] = np.where(below[:, :, p], fill[p], raw[:, :, p])
    center, scale = np.zeros(P), np.ones(P)
    for p in range(P):
        vals = raw[fitrows, :, p]
        vals = vals[~np.isnan(vals)]
        if vals.size:
            center[p] = vals.mean()
            scale[p] = vals.std() or 1.0
    X = np.where(np.isnan(raw), 0.0, (raw - center) / scale)
    y = np.array([K if r["label"] == NONE else fam_idx[safe(r["label"])] for r in keep])
    w = np.zeros(N)
    if fc.weight_by_group:
        sizes = Counter((r["kind"], r["group"]) for r, t in zip(keep, fitrows, strict=True) if t)
        for i, r in enumerate(keep):
            if fitrows[i]:
                w[i] = 1.0 / sizes[(r["kind"], r["group"])]
    else:
        w[fitrows] = 1.0
    if fc.balance_negative_sources:
        srcs = [r["source"] for r in keep]
        neg = np.array([r["kind"] == "negative" for r in keep]) & fitrows
        shuf = neg & np.array([sr == "shuffled" for sr in srcs])
        real = neg & ~shuf
        if shuf.any() and real.any():
            target = (w[shuf].sum() + w[real].sum()) / 2
            w[shuf] *= target / w[shuf].sum()
            w[real] *= target / w[real].sum()
    return Design([r["id"] for r in keep], [r["role"] for r in keep], [r["kind"] for r in keep],
                  [r["source"] for r in keep], [r["label"] for r in keep], fams, feats, raw, X, M, y, w,
                  center, scale, available, fill, fc.hmm_floor)


def export_model(d: Design, alpha: np.ndarray, B: np.ndarray, **extra) -> dict:
    """Model JSON in the form core.apply_model reads."""
    tr = d.w > 0
    counts = Counter(int(v) for v in d.y[tr])
    return {"families": d.families, "classes": d.families + [NONE], "features": d.feature_names(),
            "center": d.center.tolist(), "scale": d.scale.tolist(), "alpha": np.asarray(alpha).tolist(),
            "B": np.asarray(B).tolist(), "hmm_floor": d.hmm_floor,
            "fill": [None if np.isnan(v) else float(v) for v in (d.fill if d.fill is not None else [np.nan] * len(d.feats))],
            "available": d.available.tolist() if d.available is not None else [[True] * len(d.feats)] * len(d.families),
            "training_class_counts": {(d.families + [NONE])[k]: v for k, v in sorted(counts.items())},
            **extra}


# --------------------------------------------------------------- baseline ----

def fit_baseline(cfg: Config, d: Design) -> tuple[np.ndarray, list[dict]]:
    """HMM cutoffs from training negatives; returns cutoffs[K, P] (nan if none).

    fit.baseline_fpr is the target false-positive rate per *sequence*: with
    n profiles each cutoff uses fpr / n (Bonferroni), otherwise per-profile
    errors add up (1% on each of 50 profiles passes ~40% of negatives).
    """
    n_profiles = sum(1 for _, m in d.feats if m == "hmm") * len(d.families)
    q = 1.0 - cfg.fit.baseline_fpr / max(1, n_profiles)
    neg = (np.array(d.roles) == "train") & (d.y == len(d.families))
    cut = np.full((len(d.families), len(d.feats)), np.nan)
    rows = []
    for k, f in enumerate(d.families):
        for p, (rg, m) in enumerate(d.feats):
            if m != "hmm":
                continue
            vals = d.raw[neg, k, p]
            vals = vals[~np.isnan(vals)]
            if vals.size == 0:
                continue
            cut[k, p] = max(float(np.quantile(vals, q)), cfg.score.hmm_min_bits)
            rows.append({"family": f, "region": rg, "method": m, "cutoff_bits": round(cut[k, p], 2),
                         "negatives": int(vals.size), "quantile": round(q, 6)})
    return cut, rows


def predict_baseline(d: Design, cut: np.ndarray) -> np.ndarray:
    excess = np.where(np.isnan(cut)[None] | np.isnan(d.raw), 0.0, np.maximum(d.raw - np.nan_to_num(cut)[None], 0.0))
    s = excess.sum(axis=2)                      # N x K
    pred = s.argmax(axis=1)
    pred[s.max(axis=1) <= 0] = len(d.families)  # NONE
    return pred


# -------------------------------------------------------------------- MAP ----

def _unpack(theta: np.ndarray, K: int, P: int):
    mu = theta[0]
    alpha = theta[1:1 + K]
    beta = theta[1 + K:1 + K + P]
    delta = theta[1 + K + P:1 + K + P + K * P].reshape(K, P)
    gamma = theta[1 + K + P + K * P:]
    return mu, alpha, beta, delta, gamma


def _eta(X, M, alpha, B, gamma):
    eta = alpha[None, :] + np.einsum("nkp,kp->nk", X, B) + np.einsum("nkp,p->nk", M, gamma)
    return np.concatenate([eta, np.zeros((eta.shape[0], 1))], axis=1)


def map_objective(theta, X, M, y, w, K, P, sb, tau, sa):
    tau = np.broadcast_to(np.asarray(tau, dtype=float), (P,))[None, :]
    mu, alpha, beta, delta, gamma = _unpack(theta, K, P)
    B = beta[None, :] + delta
    full = _eta(X, M, alpha, B, gamma)
    lse = logsumexp(full, axis=1)
    nll = -np.sum(w * (full[np.arange(len(y)), y] - lse))
    pen = (np.sum(beta ** 2) / (2 * sb ** 2) + np.sum(delta ** 2 / (2 * tau ** 2))
           + np.sum((alpha - mu) ** 2) / (2 * sa ** 2) + mu ** 2 / (2 * sa ** 2) + np.sum(gamma ** 2) / 2)
    prob = np.exp(full - lse[:, None])
    onehot = np.zeros_like(prob)
    onehot[np.arange(len(y)), y] = 1.0
    R = w[:, None] * (prob - onehot)[:, :K]
    dB = np.einsum("nk,nkp->kp", R, X)
    g_alpha = R.sum(0) + (alpha - mu) / sa ** 2
    g_mu = -np.sum(alpha - mu) / sa ** 2 + mu / sa ** 2
    g_beta = dB.sum(0) + beta / sb ** 2
    g_delta = dB + delta / tau ** 2
    g_gamma = np.einsum("nk,nkp->p", R, M) + gamma
    grad = np.concatenate([[g_mu], g_alpha, g_beta, g_delta.ravel(), g_gamma])
    return nll + pen, grad


def fit_map(cfg: Config, d: Design, tau=None) -> dict:
    """Penalised fit; tau: scalar or one value per feature (default fit.tau_fixed)."""
    fc = cfg.fit
    tau = fc.tau_fixed if tau is None else tau
    tr = d.w > 0
    K, P = len(d.families), len(d.feats)
    theta0 = np.zeros(1 + K + P + K * P + P)
    res = minimize(map_objective, theta0, jac=True, method="L-BFGS-B",
                   args=(d.X[tr], d.M[tr], d.y[tr], d.w[tr], K, P, fc.beta_prior_sd, tau, fc.alpha_prior_sd),
                   options={"maxiter": 5000, "gtol": 1e-6})
    if not res.success:
        log.warning("MAP optimisation: %s", res.message)
    mu, alpha, beta, delta, gamma = _unpack(res.x, K, P)
    return {"alpha": alpha, "beta": beta, "B": beta[None, :] + delta, "gamma": gamma, "mu_alpha": mu,
            "converged": bool(res.success), "iterations": int(res.nit), "objective": float(res.fun)}


def predict_proba(d: Design, alpha: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Class probabilities (N x K+1); missing-profile cells play no role at prediction."""
    return softmax(_eta(d.X, np.zeros_like(d.M), alpha, B, np.zeros(B.shape[1])), axis=1)


# ------------------------------------------------------------------ Bayes ----

def cmdstan_available(cfg: Config) -> bool:
    try:
        import cmdstanpy
        if cfg.fit.cmdstan:
            cmdstanpy.set_cmdstan_path(cfg.fit.cmdstan)
        cmdstanpy.cmdstan_path()
        return True
    except Exception:
        return False


def plan_cores(total: int, chains: int) -> tuple[int, int]:
    """Spread `total` cores over chains: (parallel_chains, threads_per_chain).

    Chains run side by side first; leftover cores become threads inside each
    chain (reduce_sum over sequences). 16 cores, 4 chains -> 4 x 4; 2 cores,
    4 chains -> 2 chains at a time x 1 thread.
    """
    total = max(1, total)
    parallel = max(1, min(chains, total))
    return parallel, max(1, total // parallel)


def fit_bayes(cfg: Config, d: Design, workdir: Path, cache: Path) -> dict:
    import os
    import time

    import cmdstanpy

    fc = cfg.fit
    if fc.cmdstan:
        cmdstanpy.set_cmdstan_path(fc.cmdstan)
    # threaded build in its own folder, so an older single-threaded binary is never reused
    cache = cache / "threads"
    cache.mkdir(parents=True, exist_ok=True)
    stan = cache / STAN_FILE.name
    if not stan.exists() or stan.read_text() != STAN_FILE.read_text():
        shutil.copyfile(STAN_FILE, stan)
    model = cmdstanpy.CmdStanModel(stan_file=str(stan), cpp_options={"STAN_THREADS": True})
    total = fc.cores or os.cpu_count() or 1
    parallel_chains, threads = plan_cores(total, fc.chains)
    tr = d.w > 0
    n = int(tr.sum())
    grainsize = fc.grainsize or max(1, n // (4 * threads))
    data = {"N": n, "K": len(d.families), "P": len(d.feats), "y": (d.y[tr] + 1).tolist(),
            "X": d.X[tr].transpose(1, 0, 2).tolist(), "M": d.M[tr].transpose(1, 0, 2).tolist(), "w": d.w[tr].tolist(),
            "beta_prior_sd": fc.beta_prior_sd, "tau_prior_sd": fc.tau_prior_sd, "alpha_prior_sd": fc.alpha_prior_sd,
            "grainsize": grainsize}
    workdir.mkdir(parents=True, exist_ok=True)
    log.info("  Stan: %d chains, %d at a time x %d threads each (%d cores), grainsize %d",
             fc.chains, parallel_chains, threads, total, grainsize)
    t0 = time.time()
    fit = model.sample(data=data, chains=fc.chains, parallel_chains=parallel_chains, threads_per_chain=threads,
                       iter_warmup=fc.iter_warmup, iter_sampling=fc.iter_sampling, seed=fc.seed,
                       adapt_delta=fc.adapt_delta, show_progress=False, output_dir=str(workdir))
    seconds = time.time() - t0
    summ = fit.summary()
    rhat_col = next(c for c in summ.columns if c.lower().startswith("r_hat") or c.lower() == "rhat")
    ess_col = next(c for c in summ.columns if "ess_bulk" in c.lower() or c.lower() == "n_eff")
    params = summ.loc[[i for i in summ.index if not i.startswith("lp__")]]
    diag = fit.method_variables()
    divergent = int(np.sum(diag["divergent__"])) if "divergent__" in diag else -1
    alpha = fit.stan_variable("alpha")          # draws x K
    B = fit.stan_variable("B")                  # draws x K x P
    return {"alpha_draws": alpha, "B_draws": B, "summary": summ,
            "max_rhat": float(params[rhat_col].max()), "min_ess_bulk": float(params[ess_col].min()),
            "divergences": divergent, "draws": int(alpha.shape[0]),
            "beta_mean": fit.stan_variable("beta").mean(0), "tau_mean": fit.stan_variable("tau").mean(0),
            "parallel_chains": parallel_chains, "threads_per_chain": threads, "seconds": round(seconds, 1)}


def predict_bayes(d: Design, alpha_draws: np.ndarray, B_draws: np.ndarray, max_draws: int = 1000) -> np.ndarray:
    idx = np.linspace(0, alpha_draws.shape[0] - 1, min(max_draws, alpha_draws.shape[0])).astype(int)
    probs = np.zeros((len(d.ids), len(d.families) + 1))
    for s in idx:
        probs += predict_proba(d, alpha_draws[s], B_draws[s])
    return probs / len(idx)


# ---------------------------------------------------------------- metrics ----

def metrics(d: Design, pred: np.ndarray, mask: np.ndarray) -> dict:
    names = d.families + [NONE]
    y, p = d.y[mask], pred[mask]
    if y.size == 0:
        return {}
    present = sorted(set(y.tolist()))
    rec, f1 = {}, {}
    for c in present:
        tp = int(np.sum((p == c) & (y == c)))
        fn = int(np.sum((p != c) & (y == c)))
        fp = int(np.sum((p == c) & (y != c)))
        r = tp / (tp + fn) if tp + fn else 0.0
        pr = tp / (tp + fp) if tp + fp else 0.0
        rec[names[c]] = r
        f1[names[c]] = 2 * pr * r / (pr + r) if pr + r else 0.0
    K = len(d.families)
    pos, neg = y < K, y == K
    return {
        "n": int(y.size), "accuracy": float(np.mean(p == y)),
        "balanced_accuracy": float(np.mean(list(rec.values()))), "macro_f1": float(np.mean(list(f1.values()))),
        "positives_called_none": float(np.mean(p[pos] == K)) if pos.any() else None,
        "positives_wrong_family": float(np.mean((p[pos] != y[pos]) & (p[pos] != K))) if pos.any() else None,
        "negatives_called_family": float(np.mean(p[neg] != K)) if neg.any() else None,
        "recall": rec, "f1": f1,
    }


# ------------------------------------------------------------- Bayes stage --

# Bump only when the sampled posterior would change for identical inputs.
BAYES_FORMAT = 2
# [fit] settings that change the posterior beyond the design itself (the
# design - features used, floor, weights - is hashed directly). Cores,
# grainsize, CmdStan path, baseline and penalised-fit settings are left out:
# changing them never triggers a new (long) sampling run.
BAYES_KEYS = ("beta_prior_sd", "tau_prior_sd", "alpha_prior_sd", "chains", "iter_warmup", "iter_sampling",
              "adapt_delta", "seed")


def design_hash(d: Design) -> str:
    """Hash of exactly what the Bayesian model sees (training rows only)."""
    import hashlib
    tr = d.w > 0
    h = hashlib.sha256()
    h.update(repr((d.families, d.feats)).encode())
    for a in (d.X[tr], d.M[tr], d.y[tr], d.w[tr]):
        h.update(np.ascontiguousarray(np.round(np.asarray(a, dtype=float), 8)).tobytes())
    return h.hexdigest()


def bayes_core_key(cfg: Config, d: Design) -> str:
    """Identity of a posterior: training design, model settings and Stan program."""
    from .io import sha256_file, sha256_obj
    sec = cfg.section("fit")
    return sha256_obj({"stage": "bayes", "format": BAYES_FORMAT, "design": design_hash(d),
                       "model": {k: sec[k] for k in BAYES_KEYS}, "stan": sha256_file(STAN_FILE)})


def bayes_wanted(cfg: Config) -> bool:
    return cfg.fit.bayes == "true" or (cfg.fit.bayes == "auto" and cmdstan_available(cfg))


def load_bayes(cfg: Config, rdir: Path, d: Design | None = None) -> dict | None:
    """Stored posterior, only if it was sampled for the current design and settings."""
    meta_p, draws_p = rdir / "fit" / "model_bayes.json", rdir / "fit" / "bayes_draws.npz"
    if not meta_p.exists() or not draws_p.exists():
        return None
    from .io import read_json
    meta = read_json(meta_p)
    try:
        d = d if d is not None else load_design(cfg, rdir)
    except NothingToFit:
        return None
    if meta.get("key") != bayes_core_key(cfg, d):
        return None
    z = np.load(draws_p)
    return {"meta": meta, "alpha": z["alpha"], "B": z["B"]}


def run_bayes(cfg: Config, rdir: Path, cache: Path, legacy_match: bool) -> dict:
    """Sample the hierarchical model, unless an identical posterior is already on disk.

    reused   model_bayes.json carries the current key: nothing to do
    adopted  posterior written by the former single 'fit' stage, whose stamp
             matches the current data and settings: kept, key recorded
    sampled  new Stan run
    skipped  Bayesian fit disabled or CmdStan not available
    """
    from .io import read_json
    out_dir = rdir / "fit"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        d = load_design(cfg, rdir)
    except NothingToFit as e:
        return {"status": "skipped", "reason": str(e)}
    key = bayes_core_key(cfg, d)
    if load_bayes(cfg, rdir, d) is not None:
        meta = read_json(out_dir / "model_bayes.json")
        return {"status": "reused", "key": key, "diagnostics": meta.get("diagnostics")}
    meta_p, draws_p = out_dir / "model_bayes.json", out_dir / "bayes_draws.npz"
    if legacy_match and meta_p.exists() and draws_p.exists():
        meta = read_json(meta_p)
        meta["key"] = key
        write_json(meta_p, meta)
        log.info("  existing posterior adopted (sampled by the previous fit stage); no new sampling")
        return {"status": "adopted", "key": key, "diagnostics": meta.get("diagnostics")}
    if not bayes_wanted(cfg):
        if cfg.fit.bayes == "true":
            raise FitError("fit.bayes = true but CmdStan was not found (install cmdstanpy + CmdStan, or set fit.cmdstan)")
        reason = "fit.bayes = false" if cfg.fit.bayes == "false" else "CmdStan not found"
        return {"status": "skipped", "key": key, "reason": reason}

    if not (d.w > 0).any():
        raise NothingToFit("no training rows")
    log.info("  sampling the hierarchical model with Stan (%d chains x %d draws)...",
             cfg.fit.chains, cfg.fit.iter_sampling)
    b = fit_bayes(cfg, d, out_dir / "stan", cache)
    b["summary"].to_csv(out_dir / "bayes_summary.csv")
    np.savez_compressed(draws_p, alpha=b["alpha_draws"], B=b["B_draws"])
    diag = {"max_rhat": b["max_rhat"], "min_ess_bulk": b["min_ess_bulk"], "divergences": b["divergences"],
            "draws": b["draws"], "seconds": b["seconds"], "parallel_chains": b["parallel_chains"],
            "threads_per_chain": b["threads_per_chain"]}
    write_json(meta_p, export_model(
        d, b["alpha_draws"].mean(0), b["B_draws"].mean(0), model="conditional logit, hierarchical (Stan)",
        key=key, beta_mean=b["beta_mean"].tolist(), tau_mean=b["tau_mean"].tolist(),
        draws_file="bayes_draws.npz", settings=cfg.section("fit"), diagnostics=diag))
    return {"status": "sampled", "key": key, "diagnostics": diag}


# --------------------------------------------------------------- fit stage --

def load_design(cfg: Config, rdir: Path, fit_roles: tuple[str, ...] = ("train",)) -> Design:
    rows = read_tsv(rdir / "features.tsv")
    with open(rdir / "features.tsv", encoding="utf-8") as fh:
        fields = fh.readline().rstrip("\n").split("\t")
    return make_design(cfg, rows, fields, fit_roles)


def posterior_tau(cfg: Config, rdir: Path, feats: list[str], d: Design | None = None) -> tuple[np.ndarray, str]:
    """Per-feature shrinkage for the penalised fit: Bayesian posterior means when
    available and fit.tau_from_posterior, else fit.tau_fixed."""
    tau = np.full(len(feats), cfg.fit.tau_fixed)
    if not cfg.fit.tau_from_posterior:
        return tau, f"tau_fixed = {cfg.fit.tau_fixed:g}"
    b = load_bayes(cfg, rdir, d)
    if b is None:
        return tau, f"tau_fixed = {cfg.fit.tau_fixed:g} (no posterior for this design)"
    post = dict(zip(b["meta"]["features"], b["meta"]["tau_mean"], strict=True))
    tau = np.array([post.get(f, cfg.fit.tau_fixed) for f in feats])
    return tau, "posterior means of tau"


def run_fit(cfg: Config, rdir: Path) -> dict:
    """Baseline and penalised fit (seconds); Bayesian predictions from stored draws."""
    d = load_design(cfg, rdir)
    K = len(d.families)
    if not (d.w > 0).any():
        raise NothingToFit("no training rows")
    log.info("  design: %d sequences (%d train), %d families + NONE, %d features per family: %s",
             len(d.ids), int((d.w > 0).sum()), K, len(d.feats), ", ".join(f"{r}:{m}" for r, m in d.feats))

    cut, cut_rows = fit_baseline(cfg, d)
    preds = {"baseline": predict_baseline(d, cut)}
    probs = {}

    tau, tau_src = posterior_tau(cfg, rdir, d.feature_names(), d)
    m = fit_map(cfg, d, tau)
    probs["map"] = predict_proba(d, m["alpha"], m["B"])
    preds["map"] = probs["map"].argmax(1)
    out_dir = rdir / "fit"
    write_json(out_dir / "model_map.json", export_model(
        d, m["alpha"], m["B"], model="conditional logit, penalised", beta=m["beta"].tolist(),
        tau=np.asarray(tau).tolist(), tau_source=tau_src, settings=cfg.section("fit"),
        converged=m["converged"], iterations=m["iterations"], fitted_on="train"))

    ship_info = None
    if cfg.fit.ship_all_data:
        d_all = load_design(cfg, rdir, ("train", "test", "temporal_test"))
        tau_all, _ = posterior_tau(cfg, rdir, d_all.feature_names(), d)
        ms = fit_map(cfg, d_all, tau_all)
        write_json(out_dir / "model_ship.json", export_model(
            d_all, ms["alpha"], ms["B"], model="conditional logit, penalised, refitted on all labelled data",
            beta=ms["beta"].tolist(), tau=np.asarray(tau_all).tolist(), tau_source=tau_src,
            settings=cfg.section("fit"), converged=ms["converged"], iterations=ms["iterations"],
            fitted_on="train+test"))
        ship_info = {"sequences": int((d_all.w > 0).sum()), "converged": ms["converged"]}

    bayes_info = None
    b = load_bayes(cfg, rdir, d)
    if b is not None:
        if b["meta"]["families"] != d.families or b["meta"]["features"] != [f"{r}:{mm}" for r, mm in d.feats]:
            raise FitError("stored posterior does not match the current design; rerun with --from bayes")
        probs["bayes"] = predict_bayes(d, b["alpha"], b["B"])
        preds["bayes"] = probs["bayes"].argmax(1)
        meta = b["meta"]
        feats = meta["features"]
        bayes_info = dict(meta.get("diagnostics") or {})
        bayes_info["beta"] = dict(zip(feats, np.round(meta["beta_mean"], 2).tolist(), strict=True))
        bayes_info["tau"] = dict(zip(feats, np.round(meta["tau_mean"], 2).tolist(), strict=True))

    names = d.families + [NONE]
    roles = np.array(d.roles)
    pred_rows = []
    for i in np.where(roles != "train")[0]:
        row = {"id": d.ids[i], "role": d.roles[i], "kind": d.kinds[i], "source": d.sources[i], "label": d.labels[i]}
        for meth, pr in preds.items():
            row[f"{meth}_pred"] = names[pr[i]]
            if meth in probs:
                pm = probs[meth][i]
                k = int(pm[:K].argmax())
                row[f"{meth}_prob"] = round(float(pm[pr[i]]), 4)
                row[f"{meth}_prob_true"] = round(float(pm[d.y[i]]), 4)
                row[f"{meth}_p_none"] = round(float(pm[K]), 4)
                row[f"{meth}_best_family"] = names[k]
                row[f"{meth}_best_family_conditional"] = round(float(pm[k] / max(pm[:K].sum(), 1e-300)), 4)
        pred_rows.append(row)
    pred_fields = ["id", "role", "kind", "source", "label"]
    for meth in preds:
        pred_fields += [f"{meth}_pred"] + ([f"{meth}_prob", f"{meth}_prob_true", f"{meth}_p_none",
                                            f"{meth}_best_family", f"{meth}_best_family_conditional"]
                                           if meth in probs else [])
    results = {}
    for split in ("test", "temporal_test"):
        mask = roles == split
        if mask.any():
            results[split] = {meth: metrics(d, pr, mask) for meth, pr in preds.items()}
            for meth, pr in preds.items():
                by = defaultdict(list)
                for i in np.where(mask & (d.y == K))[0]:
                    by[d.sources[i]].append(pr[i] != K)
                results[split][meth]["negatives_called_family_by_source"] = {s: float(np.mean(v)) for s, v in by.items()}
    return {"design": d, "cutoffs": cut_rows, "predictions": pred_rows, "pred_fields": pred_fields,
            "results": results, "map": {k: m[k] for k in ("converged", "iterations", "objective")} | {"tau": tau_src},
            "bayes": bayes_info, "ship": ship_info}


def render_fit(cfg: Config, release: str, out: dict) -> str:
    from .audit import _table

    d: Design = out["design"]
    lines = [f"# Classifier fit: UniProt {release}", "",
             f"- Training: {int((d.w > 0).sum())} sequences (case weights 1/group size), "
             f"{len(d.families)} families + NONE; features per family: "
             + ", ".join(f"{r}:{m}" for r, m in d.feats) + ".",
             f"- baseline: per-profile HMM cutoffs from training negatives, for a {cfg.fit.baseline_fpr:.0%} "
             "false-positive rate per sequence (Bonferroni over profiles); highest summed excess wins.",
             f"- map: penalised conditional logit (shrinkage: {out['map']['tau']}); "
             f"converged: {out['map']['converged']} ({out['map']['iterations']} iterations).",
             f"- HMM floor: {cfg.fit.hmm_floor:g} bits; regions: {', '.join(cfg.fit.regions)}; "
             f"feature types: {', '.join(cfg.fit.methods)}.",
             ]
    if out.get("ship"):
        lines.append(f"- shipped model (`fit/model_ship.json`): the penalised fit refitted on all "
                     f"{out['ship']['sequences']} labelled sequences (train + test), same shrinkage; "
                     f"converged: {out['ship']['converged']}. Performance reported here is from the "
                     "train-only fit on held-out data.")
    b = out["bayes"]
    if b:
        ok = b.get("max_rhat", 9) < 1.01 and b.get("divergences", 1) == 0
        timing = (f" in {b['seconds']:.0f} s ({b['parallel_chains']} chains at a time x {b['threads_per_chain']} threads)"
                  if "seconds" in b else "")
        lines += [f"- bayes: hierarchical model in Stan, {b.get('draws', '?')} draws{timing}; "
                  f"max R-hat {b.get('max_rhat', float('nan')):.3f}, "
                  f"min bulk ESS {b.get('min_ess_bulk', float('nan')):.0f}, divergences {b.get('divergences', '?')}"
                  + ("" if ok else "  **check convergence before using these results**") + ".",
                  "", "Shared weights (posterior mean, standardised scores) and between-family spread:", "",
                  _table([{"feature": f, "beta": b["beta"][f], "tau": b["tau"][f]} for f in b["beta"]],
                         ["feature", "beta", "tau"])]
    else:
        lines.append("- bayes: no posterior for the current data and settings (see fit/bayes_status.json).")
    for split, res in out["results"].items():
        lines += ["", f"## {split.replace('_', ' ').capitalize()} set", ""]
        rows = []
        for meth, m in res.items():
            if not m:
                continue
            rows.append({"method": meth, "n": m["n"], "accuracy": f"{m['accuracy']:.3f}",
                         "balanced accuracy": f"{m['balanced_accuracy']:.3f}", "macro F1": f"{m['macro_f1']:.3f}",
                         "positives missed (NONE)": _pct(m["positives_called_none"]),
                         "positives wrong family": _pct(m["positives_wrong_family"]),
                         "negatives called a family": _pct(m["negatives_called_family"]),
                         "by negative source": ", ".join(f"{s} {_pct(v)}" for s, v in sorted(m["negatives_called_family_by_source"].items()))})
        lines.append(_table(rows, list(rows[0]) if rows else ["method"]))
        meths = [mm for mm in res if res[mm]]
        fams = sorted({f for mm in meths for f in res[mm]["recall"]})
        lines += ["Recall per class:", "",
                  _table([{"class": f, **{mm: f"{res[mm]['recall'].get(f, float('nan')):.2f}" for mm in meths}} for f in fams],
                         ["class"] + meths)]
    lines += ["", "Detailed evaluation (by sequence type, identity to training, calibration, ablations): `evaluation.md`."]
    return "\n".join(lines) + "\n"


def _pct(v) -> str:
    return "-" if v is None else f"{100 * v:.1f}%"

