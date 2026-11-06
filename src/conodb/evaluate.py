"""Evaluate: how good the classifiers are, where they fail, and why.

Reads the outputs of score/fit/bayes (never refits the Bayesian model) and
reports, on the test and temporal_test sets:

  by sequence type    complete precursors / signal present but incomplete /
                      no signal peptide (mature peptide, pro+mature fragments);
                      negatives by source
  by identity         accuracy by the sequence's highest identity to any
                      training conopeptide (evaluate.identity_bands)
  calibration         reliability of the predicted probability, Brier score,
                      log loss
  family if known     for positives, the best family when NONE is excluded
                      (the question a peptidomics user asks)
  misses              every wrongly classified positive, with its scores
  ablations           test-set results of every classifier variant of the
                      grid (`variants.classifier_variants`, the same grid the
                      select stage compares by cross-validation), penalised
                      fits, with bootstrap intervals over entry groups. The
                      variant recommended by cross-validation is marked; these
                      numbers are reported after the choice, never used for it.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np

from .build import safe
from .cluster import max_identity
from .config import Config
from .fit import (
    NONE,
    Design,
    NothingToFit,
    fit_baseline,
    fit_map,
    load_bayes,
    load_design,
    make_design,
    posterior_tau,
    predict_baseline,
    predict_bayes,
    predict_proba,
)
from .io import read_json, read_tsv

log = logging.getLogger(__name__)

METRIC_FIELDS = [
    "split",
    "method",
    "group",
    "subset",
    "n",
    "accuracy",
    "positives_correct",
    "family_if_known",
    "positives_called_none",
    "positives_wrong_family",
    "negatives_called_family",
    "brier",
    "log_loss",
]
MISS_FIELDS = [
    "id",
    "split",
    "label",
    "method",
    "predicted",
    "p_true",
    "best_family",
    "best_family_share",
    "type",
    "length",
    "regions",
    "max_identity",
    "nearest_training",
    "own_SIG_hmm",
    "own_MAT_hmm",
]
CALIB_FIELDS = ["split", "method", "bin", "n", "mean_confidence", "accuracy"]
ABLATION_FIELDS = [
    "variant",
    "note",
    "recommended",
    "features",
    "hmm_floor",
    "split",
    "subset",
    "n",
    "accuracy",
    "positives_correct",
    "positives_correct_lo",
    "positives_correct_hi",
    "family_if_known",
    "negatives_called_family",
    "macro_recall",
    "log_loss",
    "log_loss_lo",
    "log_loss_hi",
]


# ------------------------------------------------------------------ inputs --


def sequence_types(entries: list[dict], regions: list[dict], d: Design) -> list[str]:
    regs = defaultdict(set)
    for r in regions:
        regs[r["accession"]].add(r["region"])
    ent = {e["accession"]: e for e in entries}
    out = []
    for i, rid in enumerate(d.ids):
        if d.kinds[i] == "negative":
            out.append(f"negative: {d.sources[i]}")
            continue
        e = ent.get(rid)
        if e and e["complete_precursor"] == "true":
            out.append("complete precursor")
        elif "SIG" in regs.get(rid, set()):
            out.append("signal, incomplete")
        else:
            out.append("no signal peptide")
    return out


def all_sequences(entries: list[dict], negatives: list[dict]) -> dict[str, str]:
    seqs = {e["accession"]: e["sequence"] for e in entries}
    seqs.update({n["id"]: n["sequence"] for n in negatives})
    return seqs


# ---------------------------------------------------------------- metrics ---


def _summary(
    d: Design, pred: np.ndarray, prob: np.ndarray | None, mask: np.ndarray
) -> dict:
    K = len(d.families)
    y, p = d.y[mask], pred[mask]
    if y.size == 0:
        return {}
    pos, neg = y < K, y == K
    out = {"n": int(y.size), "accuracy": float(np.mean(p == y))}
    if pos.any():
        out["positives_correct"] = float(np.mean(p[pos] == y[pos]))
        out["positives_called_none"] = float(np.mean(p[pos] == K))
        out["positives_wrong_family"] = float(
            np.mean((p[pos] != y[pos]) & (p[pos] != K))
        )
        if prob is not None:
            fam = prob[mask][pos][:, :K].argmax(1)
            out["family_if_known"] = float(np.mean(fam == y[pos]))
        recalls = [float(np.mean(p[y == c] == c)) for c in sorted(set(y.tolist()))]
        out["macro_recall"] = float(np.mean(recalls))
    if neg.any():
        out["negatives_called_family"] = float(np.mean(p[neg] != K))
    if prob is not None:
        pm = np.clip(prob[mask], 1e-12, 1)
        onehot = np.zeros_like(pm)
        onehot[np.arange(len(y)), y] = 1
        out["brier"] = float(np.mean(np.sum((pm - onehot) ** 2, axis=1)))
        out["log_loss"] = float(-np.mean(np.log(pm[np.arange(len(y)), y])))
    return out


def _bands(edges: list[float]) -> list[tuple[float, float, str]]:
    e = [0.0] + sorted(edges) + [1.0000001]
    return [
        (lo, hi, f"[{lo:g}, {min(hi, 1):g}{']' if hi > 1 else ')'}")
        for lo, hi in zip(e, e[1:], strict=False)
    ]


# --------------------------------------------------------------------- run --


def run_evaluate(cfg: Config, rdir: Path) -> dict:
    ec = cfg.evaluate
    d = load_design(cfg, rdir)
    K = len(d.families)
    entries = read_tsv(rdir / "entries.tsv")
    regions = read_tsv(rdir / "regions.tsv")
    negatives = read_tsv(rdir / "negatives.tsv")
    seqs = all_sequences(entries, negatives)
    roles = np.array(d.roles)
    splits = [s for s in ("test", "temporal_test") if (roles == s).any()]
    evalmask = np.isin(roles, splits)

    # probabilities of each method (no refit of the Bayesian model)
    probs: dict[str, np.ndarray | None] = {}
    preds: dict[str, np.ndarray] = {}
    cut, _ = fit_baseline(cfg, d)
    preds["baseline"], probs["baseline"] = predict_baseline(d, cut), None
    mm = read_json(rdir / "fit" / "model_map.json")
    if mm.get("families") != d.families:
        raise RuntimeError(
            "fit/model_map.json does not match features.tsv; rerun the fit stage"
        )
    if mm.get("features") != d.feature_names():
        raise RuntimeError(
            "fit/model_map.json features differ from the current settings; rerun the fit stage"
        )
    probs["map"] = predict_proba(d, np.array(mm["alpha"]), np.array(mm["B"]))
    preds["map"] = probs["map"].argmax(1)
    b = load_bayes(cfg, rdir, d)
    if b is not None:
        probs["bayes"] = predict_bayes(d, b["alpha"], b["B"])
        preds["bayes"] = probs["bayes"].argmax(1)

    types = np.array(sequence_types(entries, regions, d))

    # identity to the training conopeptides
    train_pos = [
        (d.ids[i], seqs[d.ids[i]])
        for i in range(len(d.ids))
        if d.roles[i] == "train" and d.kinds[i] == "positive"
    ]
    queries = [(d.ids[i], seqs[d.ids[i]]) for i in np.where(evalmask)[0]]
    near = max_identity(queries, train_pos, cfg.cluster.threads)
    ident = np.array(
        [
            near.get(i, ("", np.nan))[1] if evalmask[k] else np.nan
            for k, i in enumerate(d.ids)
        ]
    )
    bands = _bands(ec.identity_bands)
    band_of = np.array(
        [
            "" if np.isnan(v) else next(lbl for lo, hi, lbl in bands if lo <= v < hi)
            for v in ident
        ]
    )

    metric_rows, calib_rows = [], []
    for split in splits:
        sm = roles == split
        for meth in preds:
            rows = [("all", "all", sm)]
            rows += [
                ("sequence type", t, sm & (types == t)) for t in sorted(set(types[sm]))
            ]
            rows += [
                ("identity to training", lbl, sm & (band_of == lbl) & (d.y < K))
                for _, _, lbl in bands
            ]
            for group, subset, mask in rows:
                if not mask.any():
                    continue
                metric_rows.append(
                    {
                        "split": split,
                        "method": meth,
                        "group": group,
                        "subset": subset,
                        **_summary(d, preds[meth], probs[meth], mask),
                    }
                )
            if probs[meth] is not None:
                conf = probs[meth].max(1)
                edges = [0.0] + sorted(ec.calibration_bins) + [1.0000001]
                for lo, hi in zip(edges, edges[1:], strict=False):
                    m = sm & (conf >= lo) & (conf < hi)
                    if m.any():
                        calib_rows.append(
                            {
                                "split": split,
                                "method": meth,
                                "bin": f"[{lo:g}, {min(hi, 1):g})",
                                "n": int(m.sum()),
                                "mean_confidence": round(float(conf[m].mean()), 3),
                                "accuracy": round(
                                    float(np.mean(preds[meth][m] == d.y[m])), 3
                                ),
                            }
                        )

    # misses (positives), model-based methods
    names = d.families + [NONE]
    scores = read_tsv(rdir / "scores.tsv")
    own = {}
    for s in scores:
        if s["method"] == "hmm" and s["region"] in ("SIG", "MAT"):
            own[(s["id"], safe(s["superfamily"]), s["region"], s["library"])] = s[
                "score"
            ]
    regs = defaultdict(list)
    for r in regions:
        regs[r["accession"]].append(r["region"])
    miss_rows = []
    for meth in [m for m in preds if probs[m] is not None]:
        for i in np.where(evalmask & (d.y < K) & (preds[meth] != d.y))[0]:
            fam = probs[meth][i, :K]
            bf = int(fam.argmax())
            lab = safe(d.labels[i])
            miss_rows.append(
                {
                    "id": d.ids[i],
                    "split": d.roles[i],
                    "label": d.labels[i],
                    "method": meth,
                    "predicted": names[preds[meth][i]],
                    "p_true": round(float(probs[meth][i, d.y[i]]), 4),
                    "best_family": names[bf],
                    "best_family_share": round(float(fam[bf] / fam.sum()), 3),
                    "type": types[i],
                    "length": len(seqs.get(d.ids[i], "")),
                    "regions": ";".join(regs.get(d.ids[i], [])),
                    "max_identity": round(float(ident[i]), 3),
                    "nearest_training": near.get(d.ids[i], ("", 0))[0],
                    "own_SIG_hmm": own.get((d.ids[i], lab, "SIG", "final"), ""),
                    "own_MAT_hmm": own.get((d.ids[i], lab, "MAT", "final"), ""),
                }
            )

    # ablations: the variant grid, penalised fits, test-set numbers with intervals
    ablation_rows = []
    if ec.ablations:
        from .select import _ci, boot_ratio, bootstrap_counts
        from .variants import classifier_variants

        rows_all = read_tsv(rdir / "features.tsv")
        with open(rdir / "features.tsv", encoding="utf-8") as fh:
            fields = fh.readline().rstrip("\n").split("\t")
        have = {c.split(":")[2] for c in fields if c.count(":") == 2}
        rec_path = rdir / "selection_recommended.json"
        recommended = (
            read_json(rec_path).get("classifier_variant") if rec_path.exists() else None
        )
        group_of = {r["id"]: f"{r['kind']}:{r['group']}" for r in rows_all}
        type_of = dict(zip(d.ids, types, strict=True))
        fc = cfg.fit
        for v in classifier_variants(cfg, have):
            cv = replace(
                cfg,
                fit=replace(
                    fc,
                    methods=list(v["methods"]),
                    regions=list(v["regions"]),
                    hmm_floor=v["hmm_floor"],
                ),
            )
            try:
                dv = make_design(cv, rows_all, fields)
            except NothingToFit:
                continue
            tau, _ = posterior_tau(cfg, rdir, dv.feature_names(), d)
            m = fit_map(cv, dv, tau)
            pr = predict_proba(dv, m["alpha"], m["B"])
            pd_ = pr.argmax(1)
            vroles = np.array(dv.roles)
            vtypes = np.array([type_of.get(i, "") for i in dv.ids])
            for split in splits:
                sm = vroles == split
                for subset, mask in [
                    ("all", sm),
                    ("complete precursor", sm & (vtypes == "complete precursor")),
                    ("no signal peptide", sm & (vtypes == "no signal peptide")),
                ]:
                    if not mask.any():
                        continue
                    s_ = _summary(dv, pd_, pr, mask)
                    idx = np.where(mask)[0]
                    gl = [group_of[dv.ids[i]] for i in idx]
                    gmap = {g: n for n, g in enumerate(dict.fromkeys(gl))}
                    gidx = np.array([gmap[g] for g in gl])
                    cnt = bootstrap_counts(len(gmap), ec.bootstrap, ec.seed)
                    yy = dv.y[idx]
                    posm = (yy < len(dv.families)).astype(float)
                    ll = -np.log(np.clip(pr[idx, yy], 1e-12, 1))
                    corr = (pd_[idx] == yy).astype(float) * posm
                    lb = (
                        boot_ratio(ll[:, None], np.ones(len(idx)), gidx, cnt)[:, 0]
                        if cnt.size
                        else np.array([])
                    )
                    pb = (
                        boot_ratio(corr[:, None], posm, gidx, cnt)[:, 0]
                        if cnt.size and posm.any()
                        else np.array([])
                    )
                    llo, lhi = _ci(lb) if lb.size else (None, None)
                    plo, phi = _ci(pb) if pb.size else (None, None)
                    ablation_rows.append(
                        {
                            "variant": v["name"],
                            "note": v["note"],
                            "recommended": v["name"] == recommended,
                            "features": ", ".join(dv.feature_names()),
                            "hmm_floor": v["hmm_floor"],
                            "split": split,
                            "subset": subset,
                            "n": s_["n"],
                            "accuracy": s_["accuracy"],
                            "positives_correct": s_.get("positives_correct"),
                            "positives_correct_lo": plo,
                            "positives_correct_hi": phi,
                            "family_if_known": s_.get("family_if_known"),
                            "negatives_called_family": s_.get(
                                "negatives_called_family"
                            ),
                            "macro_recall": s_.get("macro_recall"),
                            "log_loss": s_.get("log_loss"),
                            "log_loss_lo": llo,
                            "log_loss_hi": lhi,
                        }
                    )

    suggestions = []
    if b is not None:
        tau = dict(zip(b["meta"]["features"], b["meta"]["tau_mean"], strict=True))
        src = (
            "the penalised fit uses these posterior means"
            if cfg.fit.tau_from_posterior
            else f"the penalised fit uses tau_fixed = {cfg.fit.tau_fixed:g}"
        )
        suggestions.append(
            "Posterior between-family spread (tau) per feature: "
            + ", ".join(f"{k} {v:.2f}" for k, v in tau.items())
            + f"; {src}."
        )
    else:
        suggestions.append(
            f"No Bayesian posterior for this design: the penalised fits use tau_fixed = "
            f"{cfg.fit.tau_fixed:g}."
        )
    return {
        "metrics": metric_rows,
        "calibration": calib_rows,
        "misses": miss_rows,
        "ablations": ablation_rows,
        "methods": list(preds),
        "splits": splits,
        "suggestions": suggestions,
    }


# ------------------------------------------------------------------ report --


def _f(v, pct=False) -> str:
    if v is None or v == "":
        return "-"
    return f"{100 * v:.1f}%" if pct else f"{v:.3f}"


def render_evaluation(cfg: Config, release: str, out: dict) -> str:
    from .audit import _table

    lines = [
        f"# Evaluation: UniProt {release}",
        "",
        f"Settings: feature types {', '.join(cfg.fit.methods)}; regions {', '.join(cfg.fit.regions)}; "
        f"HMM floor {cfg.fit.hmm_floor:g} bits. Region prediction: `anatomy.md`.",
        "",
        "All numbers are on held-out sequences (test set, and temporal test set if enabled); "
        "no model saw them. `positives correct` = right family; `family if known` = right family "
        "when NONE is excluded (a sequence known to be a conopeptide).",
        "",
    ]
    for split in out["splits"]:
        lines += [f"## {split.replace('_', ' ').capitalize()} set", ""]
        for group in ("all", "sequence type", "identity to training"):
            rows = [
                r for r in out["metrics"] if r["split"] == split and r["group"] == group
            ]
            if not rows:
                continue
            title = {
                "all": "Overall",
                "sequence type": "By sequence type",
                "identity to training": "Positives by highest identity to a training conopeptide",
            }[group]
            lines += [f"### {title}", ""]
            lines.append(
                _table(
                    [
                        {
                            "method": r["method"],
                            "subset": r["subset"],
                            "n": r["n"],
                            "accuracy": _f(r.get("accuracy")),
                            "positives correct": _f(r.get("positives_correct"), True),
                            "family if known": _f(r.get("family_if_known"), True),
                            "called NONE": _f(r.get("positives_called_none"), True),
                            "negatives called a family": _f(
                                r.get("negatives_called_family"), True
                            ),
                            "Brier": _f(r.get("brier")),
                            "log loss": _f(r.get("log_loss")),
                        }
                        for r in rows
                    ],
                    [
                        "method",
                        "subset",
                        "n",
                        "accuracy",
                        "positives correct",
                        "family if known",
                        "called NONE",
                        "negatives called a family",
                        "Brier",
                        "log loss",
                    ],
                )
            )
        cal = [r for r in out["calibration"] if r["split"] == split]
        if cal:
            lines += [
                "### Calibration",
                "",
                "Within each confidence bin, accuracy should match the mean confidence.",
                "",
                _table(cal, ["method", "bin", "n", "mean_confidence", "accuracy"]),
            ]
        ab = [r for r in out["ablations"] if r["split"] == split]
        if ab:

            def _iv(r, k, pct):
                lo, hi = r.get(f"{k}_lo"), r.get(f"{k}_hi")
                base = _f(r.get(k), pct)
                return (
                    base
                    if lo is None or hi is None
                    else f"{base} [{_f(lo, pct)}, {_f(hi, pct)}]"
                )

            lines += [
                "### Variant grid on the test set (penalised fits)",
                "",
                "Reported after the defaults were chosen by cross-validation (`selection.md`); "
                "**bold** = the recommended variant. Intervals: bootstrap over entry groups. "
                "`features` per variant in `evaluation_ablations.tsv`.",
                "",
            ]
            for subset in dict.fromkeys(r["subset"] for r in ab):
                lines += [
                    f"#### {subset}",
                    "",
                    _table(
                        [
                            {
                                "variant": (
                                    f"**{r['variant']}**"
                                    if r["recommended"]
                                    else r["variant"]
                                ),
                                "note": r["note"],
                                "n": r["n"],
                                "positives correct [95% CI]": _iv(
                                    r, "positives_correct", True
                                ),
                                "family if known": _f(r["family_if_known"], True),
                                "negatives called a family": _f(
                                    r["negatives_called_family"], True
                                ),
                                "log loss [95% CI]": _iv(r, "log_loss", False),
                            }
                            for r in ab
                            if r["subset"] == subset
                        ],
                        [
                            "variant",
                            "note",
                            "n",
                            "positives correct [95% CI]",
                            "family if known",
                            "negatives called a family",
                            "log loss [95% CI]",
                        ],
                    ),
                ]
    miss = out["misses"]
    lines += [
        f"## Misclassified positives ({len(miss)} rows; all in `evaluation_misses.tsv`)",
        "",
        _table(
            miss[:40],
            [
                "method",
                "id",
                "label",
                "predicted",
                "p_true",
                "best_family",
                "best_family_share",
                "type",
                "length",
                "max_identity",
                "own_SIG_hmm",
                "own_MAT_hmm",
            ],
        ),
    ]
    if out["suggestions"]:
        lines += ["## Notes", ""] + [f"- {s}" for s in out["suggestions"]]
    return "\n".join(lines) + "\n"
