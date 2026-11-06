"""Anatomy: benchmark of precursor region prediction on held-out sequences.

Every held-out conopeptide (test, temporal_test) is cut into SIG/PRO/MAT/POST
by every region variant of `variants.region_variants` (envelope; label
transfer with extension within/cross, overhang 0/5, and any subset of the
three lightmotif cleavage-site refinements), all with profiles, models and
site matrices built from training entries only, and with two choices of
family model: the true family (region prediction alone) and the family
predicted by the penalised classifier (best family when NONE is excluded:
the end-to-end situation). POST is predicted and scored although it is not a
classification feature.

Metrics, against UniProt annotations (coordinates only where the UniProt
boundary is EXACT):
  per boundary   signal end, mature start, mature end: exact rate; the
                 pooled `boundary_exact` (all three) is the selection metric
  per region     detection, false region, start/end exact, within 1/2/5,
                 mean absolute error
  per sequence   exact mature peptide, per-residue label accuracy
  errors         residues between predicted and true boundaries, counted
                 (e.g. a trailing G left in the mature peptide)
reported overall, by sequence type and by family.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .build import safe
from .config import Config
from .core import REGIONS, best_domains, predict_segments
from .io import read_tsv
from .score import load_full_library, load_sites
from .variants import region_kwargs, region_variants

log = logging.getLogger(__name__)

BOUNDARIES = ("signal_end", "mature_start", "mature_end")
PRED_FIELDS = [
    "id",
    "family_source",
    "family",
    "true_family",
    "variant",
    "type",
    "length",
] + [
    f"{w}_{r}_{e}" for r in REGIONS for w in ("true", "pred") for e in ("start", "end")
]
METRIC_FIELDS = [
    "family_source",
    "variant",
    "subset",
    "region",
    "n_true",
    "detected",
    "false_region",
    "n_coord",
    "start_exact",
    "start_within1",
    "start_within2",
    "start_within5",
    "start_mae",
    "end_exact",
    "end_within1",
    "end_within2",
    "end_within5",
    "end_mae",
]
SEQ_FIELDS = [
    "family_source",
    "variant",
    "subset",
    "n",
    "boundary_exact",
    "boundary_n",
    "signal_end_exact",
    "mature_start_exact",
    "mature_end_exact",
    "mature_exact",
    "mature_exact_n",
    "residue_accuracy",
]
ERROR_FIELDS = ["family_source", "variant", "boundary", "direction", "residues", "n"]


# ------------------------------------------------------------ shared scoring --


def sequence_type(entry: dict, truth: dict) -> str:
    if entry["complete_precursor"] == "true":
        return "complete precursor"
    return "signal, incomplete" if "SIG" in truth else "no signal peptide"


def truth_of(regions: list[dict]) -> tuple[dict, dict]:
    """First region of each type per accession, and whether its boundaries are EXACT."""
    truth: dict[str, dict] = defaultdict(dict)
    exact: dict[tuple[str, str], bool] = {}
    for r in regions:
        if r["region"] in REGIONS and r["index"] == "1":
            truth[r["accession"]][r["region"]] = (int(r["start"]), int(r["end"]))
            exact[(r["accession"], r["region"])] = r["exact"] == "true"
    return truth, exact


def _labels(seg: dict, n: int) -> np.ndarray:
    out = np.array(["-"] * n, dtype=object)
    for r, (s, e) in seg.items():
        out[max(1, s) - 1 : min(n, e)] = r
    return out


def score_sequence(acc: str, seq: str, seg: dict, truth: dict, exact: dict) -> dict:
    """Per-sequence outcome of one prediction: boundaries (None = not scorable),
    exact mature peptide, residue accuracy, and boundary errors."""
    n = len(seq)
    out = {"boundaries": {}, "mature_exact": None, "errors": []}
    if "SIG" in truth and exact.get((acc, "SIG")):
        t = truth["SIG"][1]
        p = seg.get("SIG", (None, None))[1]
        out["boundaries"]["signal_end"] = p == t
        if p is not None and p != t:
            out["errors"].append(("signal_end", p, t))
    if "MAT" in truth and exact.get((acc, "MAT")):
        ts, te = truth["MAT"]
        ps, pe = seg.get("MAT", (None, None))
        out["boundaries"]["mature_start"] = ps == ts
        out["boundaries"]["mature_end"] = pe == te
        out["mature_exact"] = (ps, pe) == (ts, te)
        if ps is not None and ps != ts:
            out["errors"].append(("mature_start", ps - 1, ts - 1))  # as cut positions
        if pe is not None and pe != te:
            out["errors"].append(("mature_end", pe, te))
    out["residue_accuracy"] = (
        float(np.mean(_labels(seg, n) == _labels(truth, n))) if n else 0.0
    )
    return out


def error_residues(
    seq: str, kind: str, pred_cut: int, true_cut: int
) -> tuple[str, str]:
    """Residues between a predicted and a true cut, seen from the true boundary.
    'extra' = the prediction keeps residues the true region does not have
    (for mature_end: residues after the true end kept in the mature peptide)."""
    lo, hi = sorted((pred_cut, true_cut))
    res = seq[lo:hi]
    if kind == "mature_end":
        return ("extra" if pred_cut > true_cut else "missing"), res
    if kind == "mature_start":
        return ("extra" if pred_cut < true_cut else "missing"), res
    return ("longer" if pred_cut > true_cut else "shorter"), res


def aggregate(items: list[dict]) -> dict:
    """Pooled rates over per-sequence outcomes."""
    out = {"n": len(items)}
    tot = hit = 0
    for b in BOUNDARIES:
        vals = [it["boundaries"][b] for it in items if b in it["boundaries"]]
        out[f"{b}_exact"] = _rate(sum(vals), len(vals))
        tot += len(vals)
        hit += sum(vals)
    out["boundary_exact"], out["boundary_n"] = _rate(hit, tot), tot
    mat = [it["mature_exact"] for it in items if it["mature_exact"] is not None]
    out["mature_exact"], out["mature_exact_n"] = _rate(sum(mat), len(mat)), len(mat)
    out["residue_accuracy"] = (
        round(float(np.mean([it["residue_accuracy"] for it in items])), 4)
        if items
        else ""
    )
    return out


def _rate(k: int, n: int):
    return round(k / n, 4) if n else ""


# --------------------------------------------------------------------- run ---


def run_anatomy(cfg: Config, rdir: Path) -> dict:
    an = cfg.anatomy
    entries = {e["accession"]: e for e in read_tsv(rdir / "entries.tsv")}
    part = read_tsv(rdir / "partition.tsv")
    truth, exact = truth_of(read_tsv(rdir / "regions.tsv"))
    preds = (
        {p["id"]: p for p in read_tsv(rdir / "predictions.tsv")}
        if (rdir / "predictions.tsv").exists()
        else {}
    )
    variants = region_variants(cfg)

    held = [
        p
        for p in part
        if p["scope"] == "conopeptide"
        and p["role"] in ("test", "temporal_test")
        and p["tier"] in ("1", "2")
    ]
    full_hmms, full_labels = load_full_library(rdir / "profiles" / "final")
    families = sorted({n.split("_", 1)[1].rsplit("_", 1)[0] for n in full_labels})
    from pyhmmer import plan7

    region_hmms = []
    if (rdir / "profiles" / "final.hmm").exists() and (
        rdir / "profiles" / "final.hmm"
    ).stat().st_size:
        with plan7.HMMFile(rdir / "profiles" / "final.hmm") as fh:
            region_hmms = list(fh)
    sites = load_sites(rdir / "profiles" / "final.sites.json", an.pseudocount)
    seqs = [(p["accession"], entries[p["accession"]]["sequence"]) for p in held]
    sc = cfg.score
    hits_full = best_domains(
        full_hmms,
        seqs,
        min_bits=sc.hmm_min_bits,
        max_mode=sc.max_mode,
        cpus=sc.threads,
        alignments=True,
    )
    hits_reg = best_domains(
        region_hmms,
        seqs,
        min_bits=sc.hmm_min_bits,
        max_mode=sc.max_mode,
        cpus=sc.threads,
    )

    pred_rows, items = [], []
    for p in held:
        acc = p["accession"]
        seq = entries[acc]["sequence"]
        typ = sequence_type(entries[acc], truth[acc])
        true_fam = safe(p["superfamily"])
        sources = [("true", true_fam)]
        pr = preds.get(acc)
        if pr and pr.get("map_best_family"):
            sources.append(("predicted", pr["map_best_family"]))
        for src, fam in sources:
            if fam not in families:
                continue
            env = {
                rg: (h.start, h.end)
                for rg in REGIONS
                if (h := hits_reg.get((acc, f"CONOPEP_{fam}_{rg}"))) is not None
                and h.score > 0
            }
            for v in variants:
                seg = predict_segments(
                    seq,
                    full_hit=hits_full.get((acc, f"CONOPEP_{fam}_FULL")),
                    labels=full_labels.get(f"CONOPEP_{fam}_FULL"),
                    region_hits=env,
                    sites=sites if v["sites"] else None,
                    **region_kwargs(cfg, v),
                )
                row = {
                    "id": acc,
                    "family_source": src,
                    "family": fam,
                    "true_family": true_fam,
                    "variant": v["name"],
                    "type": typ,
                    "length": len(seq),
                }
                for rg in REGIONS:
                    t, q = truth[acc].get(rg), seg.get(rg)
                    row[f"true_{rg}_start"], row[f"true_{rg}_end"] = (
                        t if t else ("", "")
                    )
                    row[f"pred_{rg}_start"], row[f"pred_{rg}_end"] = (
                        q if q else ("", "")
                    )
                pred_rows.append(row)
                it = score_sequence(acc, seq, seg, truth[acc], exact)
                it.update(
                    {
                        "id": acc,
                        "seq": seq,
                        "seg": seg,
                        "truth": truth[acc],
                        "family_source": src,
                        "variant": v["name"],
                        "type": typ,
                        "true_family": true_fam,
                    }
                )
                items.append(it)

    metric_rows, seq_rows, error_rows = [], [], []
    by = defaultdict(list)
    for it in items:
        for subset in ("all", f"type: {it['type']}", f"family: {it['true_family']}"):
            by[(it["family_source"], it["variant"], subset)].append(it)
    order = {v["name"]: i for i, v in enumerate(variants)}
    for (src, var, subset), its in sorted(
        by.items(), key=lambda kv: (kv[0][0], order[kv[0][1]], kv[0][2])
    ):
        seq_rows.append(
            {"family_source": src, "variant": var, "subset": subset, **aggregate(its)}
        )
        for rg in REGIONS:
            has_t = [it for it in its if rg in it["truth"]]
            m = {
                "family_source": src,
                "variant": var,
                "subset": subset,
                "region": rg,
                "n_true": len(has_t),
                "detected": _rate(sum(rg in it["seg"] for it in has_t), len(has_t)),
                "false_region": _rate(
                    sum(rg in it["seg"] and rg not in it["truth"] for it in its),
                    sum(rg not in it["truth"] for it in its),
                ),
            }
            coord = [
                it for it in has_t if rg in it["seg"] and exact.get((it["id"], rg))
            ]
            m["n_coord"] = len(coord)
            for side, i in (("start", 0), ("end", 1)):
                err = np.array(
                    [abs(it["seg"][rg][i] - it["truth"][rg][i]) for it in coord]
                )
                if err.size:
                    m[f"{side}_exact"] = round(float(np.mean(err == 0)), 4)
                    for k in (1, 2, 5):
                        m[f"{side}_within{k}"] = round(float(np.mean(err <= k)), 4)
                    m[f"{side}_mae"] = round(float(err.mean()), 3)
            metric_rows.append(m)
        if subset == "all":
            cnt = Counter()
            for it in its:
                for kind, pc, tc in it["errors"]:
                    d, res = error_residues(it["seq"], kind, pc, tc)
                    cnt[
                        (kind, d, res if len(res) <= 6 else f"{len(res)} residues")
                    ] += 1
            for (kind, d, res), n in sorted(
                cnt.items(), key=lambda kv: (-kv[1], kv[0])
            ):
                error_rows.append(
                    {
                        "family_source": src,
                        "variant": var,
                        "boundary": kind,
                        "direction": d,
                        "residues": res,
                        "n": n,
                    }
                )
    return {
        "predictions": pred_rows,
        "metrics": metric_rows,
        "sequences": seq_rows,
        "errors": error_rows,
        "n": len(held),
        "families": families,
        "sites": bool(sites and sites.matrices),
        "variants": variants,
    }


def render_anatomy(
    cfg: Config, release: str, out: dict, recommended: str | None = None
) -> str:
    from .audit import _table
    from .variants import region_name

    an = cfg.anatomy
    current = region_name(
        an.method, an.extension, an.terminal_overhang, an.site_types, an.cterm_trim
    )
    lines = [
        f"# Region prediction benchmark: UniProt {release}",
        "",
        f"{out['n']} held-out conopeptides; profiles, full-precursor models and cleavage-site matrices "
        "built from training entries only. Coordinates are scored only where the UniProt boundary is EXACT. "
        "These are test-set numbers: the defaults are chosen by cross-validation on training entries "
        "(`selection.md`), not from this table.",
        "",
        f"- full-precursor models: match states = columns occupied >= {an.min_occupancy:g}; labels = majority region",
        "- transfer/within: unaligned ends extended only within the region of the last aligned state; "
        "transfer/cross: across regions; ohN: N residues past where the extension stops keep its label",
        f"- sites (S = signal end, MS = mature start, ME = mature end): lightmotif matrices of {an.site_left} "
        f"residues before / {an.site_right} after each cleavage; search +/- {an.search_radius} residues, "
        f"penalty {an.distance_penalty:g} per residue shifted; Cterm = a C-terminal mature peptide loses a "
        "trailing G[KR]{0,3} (amidation)"
        + ("" if out["sites"] else " (no site matrices: too few exact boundaries)"),
        f"- current settings (framework feature, conodictor): **{current}**"
        + (
            f"; recommended by cross-validation: **{recommended}**"
            if recommended
            else ""
        ),
        "",
    ]
    for src in ("true", "predicted"):
        seq = [
            r
            for r in out["sequences"]
            if r["family_source"] == src and r["subset"] == "all"
        ]
        if not seq:
            continue
        lines += [
            f"## Family model: {src} family",
            "",
            "### Boundaries and mature peptide",
            "",
            _table(
                seq,
                [
                    "variant",
                    "n",
                    "boundary_exact",
                    "signal_end_exact",
                    "mature_start_exact",
                    "mature_end_exact",
                    "mature_exact",
                    "residue_accuracy",
                ],
            ),
        ]
        typ = [
            r
            for r in out["sequences"]
            if r["family_source"] == src and r["subset"].startswith("type:")
        ]
        if typ:
            piv = defaultdict(dict)
            for r in typ:
                piv[r["variant"]][r["subset"][6:]] = r["mature_exact"]
            cols = sorted({r["subset"][6:] for r in typ})
            lines += [
                "### Exact mature peptide by sequence type",
                "",
                _table(
                    [
                        {"variant": v, **piv[v]}
                        for v in dict.fromkeys(r["variant"] for r in typ)
                    ],
                    ["variant"] + cols,
                ),
            ]
        reg = [
            r
            for r in out["metrics"]
            if r["family_source"] == src
            and r["subset"] == "all"
            and r["variant"] in (current, recommended)
        ]
        if reg:
            lines += [
                "### Per region (current and recommended settings)",
                "",
                _table(
                    reg,
                    [
                        "variant",
                        "region",
                        "n_true",
                        "detected",
                        "false_region",
                        "n_coord",
                        "start_exact",
                        "start_within2",
                        "start_mae",
                        "end_exact",
                        "end_within2",
                        "end_mae",
                    ],
                ),
            ]
    err = [
        r
        for r in out["errors"]
        if r["family_source"] == "true" and r["variant"] == current
    ]
    if err:
        lines += [
            f"## Boundary errors of the current settings (true family; top 15 of {len(err)})",
            "",
            "`residues` lie between the predicted and the true boundary; `extra` = kept in the predicted "
            "mature peptide although outside the true one.",
            "",
            _table(err[:15], ["boundary", "direction", "residues", "n"]),
        ]
    lines += [
        "All rows, including per-family results and within-1/5 rates: `anatomy_metrics.tsv`, "
        "`anatomy_sequences.tsv`; boundary errors of every variant: `anatomy_errors.tsv`; every prediction: "
        "`anatomy_predictions.tsv`."
    ]
    return "\n".join(lines) + "\n"
