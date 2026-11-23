#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Anicet Ebou and Dominique Koua
# Dual-licensed: see LICENSING.md for commercial licensing.
"""Benchmark of conopeptide classifiers on a conodb release's test set.

    python benchmark/benchmark.py export RELEASE OUT        # test FASTA, ConoPrec chunks, truth table
    python benchmark/benchmark.py nn RELEASE OUT            # nearest-neighbour baseline
    python benchmark/benchmark.py score RELEASE OUT [--conodictor2 summary.csv] [--conosorter DIR]
                                                    [--conoprec FILE ...] [--table NAME=FILE ...]

ConoDictor 3 is scored from the release's own predictions (fit/model_map.json
on features from the training libraries: what scripts/validate_release.py in
conodictor reproduces exactly), never from the shipped database, whose
profiles and model include the test sequences.

Other tools' labels are mapped to the release's families with
benchmark/labels.tsv (tool, tool_label, family; family OTHER = a superfamily
the release does not have, NONE = not a conopeptide). Every unmapped label is
reported. Metrics, per tool and subset, with 95% bootstrap intervals over
entry groups: positives correct, positives called NONE, positives given a
wrong family, negatives called a family. Subsets: all, by sequence type,
positives of families the tool covers, and test entries first public after
the tool's release (it cannot have seen them).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOOL_DATES = {"ConoSorter": "2013-10-01", "ConoDictor 2": "2021-07-07"}   # release of each tool (name prefix)
NONE, OTHER = "NONE", "OTHER"


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def safe(name: str) -> str:                     # as conodb.build.safe
    return re.sub(r"[^A-Za-z0-9.+-]+", "-", name).strip("-") or "unnamed"


# ------------------------------------------------------------------ export ---

def truth_table(rdir: Path) -> list[dict]:
    preds = read_tsv(rdir / "predictions.tsv")
    test_ids = {p["id"] for p in preds if p["role"] in ("test", "temporal_test")}
    ent = {e["accession"]: e for e in read_tsv(rdir / "entries.tsv")}
    part = {p["accession"]: p for p in read_tsv(rdir / "partition.tsv")}
    neg = {n["id"]: n for n in read_tsv(rdir / "negatives.tsv")}
    regs = defaultdict(set)
    for r in read_tsv(rdir / "regions.tsv"):
        regs[r["accession"]].add(r["region"])
    rows, dropped = [], []
    for i in sorted(test_ids):
        if i in part and i in ent:
            e, p = ent[i], part[i]
            typ = ("complete precursor" if e["complete_precursor"] == "true"
                   else "signal, incomplete" if "SIG" in regs[i] else "no signal peptide")
            rows.append({"id": i, "kind": "positive", "label": safe(p["superfamily"]), "type": typ,
                         "group": p["group"], "first_public": e.get("first_public", ""), "length": len(e["sequence"]),
                         "sequence": e["sequence"]})
        elif i in neg:
            n = neg[i]
            rows.append({"id": i, "kind": "negative", "label": NONE, "type": f"negative: {n['source']}",
                         "group": n["group"], "first_public": "", "length": len(n["sequence"]),
                         "sequence": n["sequence"]})
        else:
            dropped.append(i)
    if dropped:
        print(f"WARNING: {len(dropped)} test sequence(s) of predictions.tsv not found in partition/entries/negatives "
              f"and left out: {', '.join(dropped[:10])}", file=sys.stderr)
    return rows


def cmd_export(rdir: Path, out: Path, chunk: int = 450) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rows = truth_table(rdir)
    write_tsv(out / "truth.tsv", rows, ["id", "kind", "label", "type", "group", "first_public", "length"])
    fa = "".join(f">{r['id']}\n{r['sequence']}\n" for r in rows)
    (out / "test.faa").write_text(fa)
    for k in range(0, len(rows), chunk):
        part = rows[k:k + chunk]
        (out / f"conoprec_{k // chunk + 1:02d}.faa").write_text("".join(f">{r['id']}\n{r['sequence']}\n" for r in part))
    n_pos = sum(r["kind"] == "positive" for r in rows)
    print(f"{len(rows)} test sequences ({n_pos} conopeptides, {len(rows) - n_pos} negatives) -> {out}/test.faa; "
          f"{(len(rows) + chunk - 1) // chunk} ConoPrec files of <= {chunk}")


# ---------------------------------------------------------------------- nn ---

def cmd_nn(rdir: Path, out: Path, threads: int = 0) -> None:
    """Best BLOSUM62 Smith-Waterman hit (pyopal) among training sequences; its label is the prediction."""
    import os
    from multiprocessing.pool import ThreadPool

    import pyopal

    out.mkdir(parents=True, exist_ok=True)
    ent = {e["accession"]: e["sequence"] for e in read_tsv(rdir / "entries.tsv")}
    neg = {n["id"]: n["sequence"] for n in read_tsv(rdir / "negatives.tsv")}
    seq = ent | neg
    with open(rdir / "features.tsv", encoding="utf-8", newline="") as fh:      # the classifier's training rows
        train = sorted((r["id"], safe(r["label"]) if r["kind"] == "positive" else NONE)
                       for r in csv.DictReader(fh, delimiter="\t") if r["role"] == "train")
    db = pyopal.Database([seq[i].upper() for i, _ in train])
    test = [r for r in truth_table(rdir)]
    rows = []
    with ThreadPool(threads or os.cpu_count() or 1) as pool:
        for r in test:
            res = pyopal.align(r["sequence"].upper(), db, mode="score", pool=pool, ordered=True)
            best = max(res, key=lambda x: (x.score, -x.target_index))
            rows.append({"id": r["id"], "prediction": train[best.target_index][1],
                         "nearest": train[best.target_index][0], "score": best.score})
    write_tsv(out / "nn_predictions.tsv", rows, ["id", "prediction", "nearest", "score"])
    print(f"nearest neighbour: {len(rows)} predictions -> {out}/nn_predictions.tsv")


# ----------------------------------------------------------------- parsers ---

def parse_conodictor2(path: Path, column: str = "definitive_pred") -> dict[str, str]:
    """summary.csv of ConoDictor 2.x (run with -a): sequence,hmm_pred,pssm_pred,definitive_pred.
    'CONFLICT X and Y' is kept as such (a conopeptide call without a single family)."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["sequence"].split()[0]] = row[column].strip()
    return out


def parse_table(path: Path) -> dict[str, str]:
    """Any two-column file: id <tab> label (header optional)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] not in ("id", "sequence"):
            out[parts[0].split()[0]] = parts[1].strip()
    return out


def _conosorter_rows(path: Path) -> list[dict]:
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not lines:
        return []
    head = lines[0].split("|")
    rows = []
    for ln in lines[1:]:
        vals = ln.split("|")
        rows.append({"id": re.sub(r"_-?\d+_$", "", vals[1]), "cols": list(zip(head, vals, strict=False))})
    return rows


def parse_conosorter(path: Path) -> dict[str, dict[str, str]]:
    """ConoSorter 1.x output folder (or prefix): <x>Regex.tab and <x>pHMM.tab.

    Superfamily of a sequence = the superfamily its signal, pro-region and
    mature hits agree on; different superfamilies -> CONFLICT; no hit or a
    sequence absent from the file -> NONE. Superfamilies ConoSorter names
    after a signal sequence (e.g. 'MTFLLLLVSV') are not in the release (OTHER).
    ConoSorter runs its pHMMs only on sequences its regular expressions left
    unclassified (every pHMM.tab row has no regex superfamily), so the tool's
    own answer is the regex call, else the pHMM call: 'ConoSorter'. The regex
    stage alone is reported too: 'ConoSorter (regex only)'."""
    path = Path(path)
    out = {}
    for kind in ("Regex", "pHMM"):
        files = sorted(path.glob(f"*{kind}.tab")) if path.is_dir() else [Path(f"{path}{kind}.tab")]
        if not files or not files[0].exists():
            raise SystemExit(f"no *{kind}.tab in {path}")
        pred = {}
        for r in _conosorter_rows(files[0]):
            fams = []
            for name, val in r["cols"]:
                if name.startswith("Superfamily (") and val not in ("-", ""):
                    fams.append(re.sub(r"\s*\((Signal|Pro-region|Mature)\)\s*$", "", val).strip())
            fams = sorted(set(fams))
            pred[r["id"]] = "" if not fams else fams[0] if len(fams) == 1 else "CONFLICT " + " and ".join(fams)
        out[kind] = pred
    combined = dict(out["Regex"])
    for i, v in out["pHMM"].items():
        if not combined.get(i):
            combined[i] = v
    return {"ConoSorter": combined, "ConoSorter (regex only)": out["Regex"]}


def parse_conoprec(paths: list[Path]) -> dict[str, str]:
    raise SystemExit("ConoPrec parser not written yet: send one result file")


# ------------------------------------------------------------------- score ---

def label_map() -> dict[tuple[str, str], str]:
    p = HERE / "labels.tsv"
    return {(r["tool"], r["tool_label"]): r["family"] for r in read_tsv(p)} if p.exists() else {}


def normalise(tool: str, raw: str, families: set[str], mapping: dict, unmapped: dict) -> str:
    raw = (raw or "").strip()
    if (tool, raw) in mapping:
        return mapping[(tool, raw)]
    if raw in ("", "UNKNOWN", "NONE", "None", "-", "NA", "unclassified"):
        return NONE
    if raw.startswith("CONFLICT"):
        return "CONFLICT"
    if raw in families:
        return raw
    if safe(raw) in families:
        return safe(raw)
    unmapped[tool][raw] += 1
    return OTHER


def covered(tool: str, families: set[str], mapping: dict, predicted: set[str]) -> set[str]:
    """Release families a tool can output (from its label map and what it predicted)."""
    from_map = {f for (t, _), f in mapping.items() if t == tool and f in families}
    return from_map | (predicted & families)


def bootstrap_diff(a: np.ndarray, b: np.ndarray, groups: np.ndarray, B: int, seed: int) -> tuple[float, float]:
    """95% interval of mean(a) - mean(b) for paired outcomes, entry groups resampled."""
    if not len(a) or B == 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    ug, gi = np.unique(groups, return_inverse=True)
    ha = np.bincount(gi, weights=a.astype(float), minlength=len(ug))
    hb = np.bincount(gi, weights=b.astype(float), minlength=len(ug))
    n = np.bincount(gi, minlength=len(ug)).astype(float)
    cnt = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=B)
    d = (cnt @ ha - cnt @ hb) / np.maximum(cnt @ n, 1)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def bootstrap_rate(hit: np.ndarray, groups: np.ndarray, B: int, rng) -> tuple[float, float]:
    if not len(hit) or B == 0:
        return (float("nan"), float("nan"))
    ug, gi = np.unique(groups, return_inverse=True)
    h = np.bincount(gi, weights=hit.astype(float), minlength=len(ug))
    n = np.bincount(gi, minlength=len(ug)).astype(float)
    cnt = rng.multinomial(len(ug), np.full(len(ug), 1 / len(ug)), size=B)
    r = (cnt @ h) / np.maximum(cnt @ n, 1)
    return float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def cmd_score(rdir: Path, out: Path, a) -> None:
    out.mkdir(parents=True, exist_ok=True)
    truth = {r["id"]: r for r in read_tsv(out / "truth.tsv")} if (out / "truth.tsv").exists() else \
        {r["id"]: r for r in truth_table(rdir)}
    preds = {p["id"]: p for p in read_tsv(rdir / "predictions.tsv")}
    families = {safe(r["label"]) for r in truth.values() if r["kind"] == "positive"} | \
        {p["map_best_family"] for p in preds.values() if p.get("map_best_family")}
    mapping = label_map()
    tools: dict[str, dict[str, str]] = {"ConoDictor 3": {i: preds[i]["map_pred"] for i in truth}}
    if any(p.get("bayes_pred") for p in preds.values()):
        tools["ConoDictor 3 (Bayes)"] = {i: preds[i]["bayes_pred"] for i in truth}
    tools["Nearest neighbour"] = {i: v for i, v in parse_table(out / "nn_predictions.tsv").items()} \
        if (out / "nn_predictions.tsv").exists() else {}
    if a.conodictor2:
        v = a.conodictor2_version
        tools[f"ConoDictor {v}"] = parse_conodictor2(a.conodictor2)
        tools[f"ConoDictor {v} (HMM only)"] = parse_conodictor2(a.conodictor2, "hmm_pred")
    if a.conosorter:
        tools |= parse_conosorter(a.conosorter)
    if a.conoprec:
        tools["ConoPrec"] = parse_conoprec(a.conoprec)
    for spec in a.table or []:
        name, path = spec.split("=", 1)
        tools[name] = parse_table(Path(path))
    tools = {k: v for k, v in tools.items() if v}

    unmapped = defaultdict(lambda: defaultdict(int))
    rng = np.random.default_rng(20261106)
    rows, per_seq = [], []
    ref = {i: normalise("ConoDictor 3", tools["ConoDictor 3"].get(i, ""), families, mapping, unmapped) for i in truth}
    for tool, raw in tools.items():
        pred = {i: normalise(tool, raw.get(i, ""), families, mapping, unmapped) for i in truth}
        cov = covered(tool, families, mapping, set(pred.values())) if not tool.startswith("ConoDictor 3") \
            and tool != "Nearest neighbour" else families
        date = next((d for k, d in TOOL_DATES.items() if tool.startswith(k)), None)
        subsets = {"all positives": lambda r: r["kind"] == "positive",
                   "complete precursor": lambda r: r["type"] == "complete precursor",
                   "signal, incomplete": lambda r: r["type"] == "signal, incomplete",
                   "no signal peptide": lambda r: r["type"] == "no signal peptide",
                   "positives, families the tool covers": lambda r, c=cov: r["kind"] == "positive" and r["label"] in c,
                   "positives first public after 2021-07-07 (unseen by every tool)":
                       lambda r: r["kind"] == "positive" and r["first_public"] > "2021-07-07",
                   "all negatives": lambda r: r["kind"] == "negative"}
        if date:
            subsets[f"positives first public after {date} (unseen by this tool)"] = \
                lambda r, d=date: r["kind"] == "positive" and r["first_public"] > d
        for t in sorted({r["type"] for r in truth.values() if r["kind"] == "negative"}):
            subsets[t] = lambda r, t=t: r["type"] == t
        for sname, f in subsets.items():
            sel = [r for r in truth.values() if f(r)]
            if not sel:
                continue
            ids = [r["id"] for r in sel]
            grp = np.array([f"{r['kind']}:{r['group']}" for r in sel])
            p = np.array([pred[i] for i in ids])
            y = np.array([r["label"] for r in sel])
            row = {"tool": tool, "subset": sname, "n": len(sel)}
            pr = np.array([ref[i] for i in ids])
            seed = int.from_bytes(sname.encode()[:8].ljust(8, b"_"), "little") % (2 ** 32)   # same resamples for every tool
            if sel[0]["kind"] == "positive":
                correct, none = p == y, p == NONE
                lo, hi = bootstrap_rate(correct, grp, a.bootstrap, rng)
                dlo, dhi = bootstrap_diff(correct, pr == y, grp, a.bootstrap, seed)
                row |= {"positives_correct": correct.mean(), "positives_correct_lo": lo, "positives_correct_hi": hi,
                        "called_none": none.mean(), "wrong_family": (~correct & ~none).mean(),
                        "diff_vs_v3": correct.mean() - (pr == y).mean(), "diff_vs_v3_lo": dlo, "diff_vs_v3_hi": dhi}
            else:
                called = p != NONE
                lo, hi = bootstrap_rate(called, grp, a.bootstrap, rng)
                dlo, dhi = bootstrap_diff(called, pr != NONE, grp, a.bootstrap, seed)
                row |= {"negatives_called_family": called.mean(), "negatives_called_family_lo": lo,
                        "negatives_called_family_hi": hi, "diff_vs_v3": called.mean() - (pr != NONE).mean(),
                        "diff_vs_v3_lo": dlo, "diff_vs_v3_hi": dhi}
            rows.append(row)
        for i in truth:
            per_seq.append({"tool": tool, "id": i, "truth": truth[i]["label"], "type": truth[i]["type"],
                            "raw": raw.get(i, ""), "prediction": pred[i]})
    fields = ["tool", "subset", "n", "positives_correct", "positives_correct_lo", "positives_correct_hi",
              "called_none", "wrong_family", "negatives_called_family", "negatives_called_family_lo",
              "negatives_called_family_hi", "diff_vs_v3", "diff_vs_v3_lo", "diff_vs_v3_hi"]
    fmt = [{k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()} for r in rows]
    write_tsv(out / "benchmark.tsv", fmt, fields)
    write_tsv(out / "benchmark_predictions.tsv", per_seq, ["tool", "id", "truth", "type", "raw", "prediction"])
    write_tsv(out / "benchmark_unmapped.tsv",
              [{"tool": t, "label": lab, "n": n} for t, d in unmapped.items() for lab, n in sorted(d.items())],
              ["tool", "label", "n"])
    (out / "benchmark.md").write_text(render(rows, rdir, tools, unmapped))
    print((out / "benchmark.md").read_text())


def _pct(v, lo=None, hi=None) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "-"
    s = f"{100 * v:.1f}%"
    if lo is not None and not np.isnan(lo):
        s += f" [{100 * lo:.1f}, {100 * hi:.1f}]"
    return s


def _dpct(r: dict) -> str:
    if r["tool"] == "ConoDictor 3":
        return "-"
    d, lo, hi = r["diff_vs_v3"], r["diff_vs_v3_lo"], r["diff_vs_v3_hi"]
    if np.isnan(lo):
        return f"{100 * d:+.1f}"
    return f"{100 * d:+.1f} [{100 * lo:+.1f}, {100 * hi:+.1f}]"


def render(rows: list[dict], rdir: Path, tools: dict, unmapped: dict) -> str:
    lines = [f"# Benchmark: test set of {rdir.name}", "",
             "ConoDictor 3 = the release's evaluation model (training libraries only). Intervals: bootstrap over "
             "entry groups; differences to ConoDictor 3 are paired (same sequences, same resamples). Other tools' labels mapped with `benchmark/labels.tsv`; CONFLICT (two superfamilies) "
             "counts as a conopeptide call without the right family; a label that is not a release superfamily "
             "(OTHER) counts as a wrong family.", ""]
    subsets = list(dict.fromkeys(r["subset"] for r in rows))
    for s in subsets:
        rs = [r for r in rows if r["subset"] == s]
        if "positives_correct" in rs[0]:
            lines += [f"## {s}", "", "| tool | n | correct [95% CI] | called NONE | wrong family | difference to ConoDictor 3 [95% CI] |",
                      "|---|---|---|---|---|---|"]
            lines += [f"| {r['tool']} | {r['n']} | {_pct(r['positives_correct'], r['positives_correct_lo'], r['positives_correct_hi'])} "
                      f"| {_pct(r['called_none'])} | {_pct(r['wrong_family'])} | {_dpct(r)} |" for r in rs]
        else:
            lines += [f"## {s}", "", "| tool | n | called a family [95% CI] | difference to ConoDictor 3 [95% CI] |",
                      "|---|---|---|---|"]
            lines += [f"| {r['tool']} | {r['n']} | {_pct(r['negatives_called_family'], r['negatives_called_family_lo'], r['negatives_called_family_hi'])} "
                      f"| {_dpct(r)} |" for r in rs]
        lines.append("")
    if any(unmapped.values()):
        lines += ["## Unmapped labels (counted as OTHER: a family the release does not have)", ""]
        lines += [f"- {t}: " + ", ".join(f"{k} ({v})" for k, v in d.items()) for t, d in unmapped.items() if d]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("export", "nn", "score"):
        s = sub.add_parser(c)
        s.add_argument("release", type=Path)
        s.add_argument("out", type=Path)
        if c == "nn":
            s.add_argument("-j", "--threads", type=int, default=0)
        if c == "score":
            s.add_argument("--conodictor2", type=Path, help="summary.csv of ConoDictor 2.x run with -a")
            s.add_argument("--conodictor2-version", default="2.4.0", help="version that wrote it [2.4.0]")
            s.add_argument("--conosorter", type=Path, help="ConoSorter output folder (Regex.tab, pHMM.tab)")
            s.add_argument("--conoprec", type=Path, nargs="+", help="ConoPrec result files")
            s.add_argument("--table", nargs="+", metavar="NAME=FILE", help="other tools: id<TAB>label files")
            s.add_argument("--bootstrap", type=int, default=1000)
    a = ap.parse_args()
    if a.cmd == "export":
        cmd_export(a.release, a.out)
    elif a.cmd == "nn":
        cmd_nn(a.release, a.out, a.threads)
    else:
        cmd_score(a.release, a.out, a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
