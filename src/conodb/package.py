"""Package: the database conodictor loads, and the supplementary bundle.

package/conodictor-db-<release>/
  profiles.hmm          region HMMs (ship set: train + test), GA cutoffs calibrated on negatives
  profiles.prf          region PSSMs (if pftools was used), /CUT_OFF calibrated on negatives
  profiles.pwm.json     lightmotif windows per region (if built)
  precursors.hmm        full-precursor HMMs (region prediction)
  precursors.labels.json  region label of every match state
  sites.json            cleavage-site windows (lightmotif matrices are built from them)
  framework.json        cysteine-framework counts per family (all labelled mature peptides)
  model.json            classifier (penalised fit refitted on all labelled data)
  settings.json         every setting needed to compute features exactly as in training
  manifest.json         versions, UniProt release, counts, sha256 of every file
  README.md
  supplementary/        every table behind the evaluation, with an index (README.md)
and conodictor-db-<release>.tar.gz next to it.

Cutoffs: for each profile, the (1 - cutoff_fpr / number of profiles)
quantile of all negatives' scores (false-positive rate cutoff_fpr per
sequence over the whole library), never below 0 bits for HMMs.
"""

from __future__ import annotations

import logging
import platform
import re
import shutil
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from . import __version__
from .build import safe
from .config import Config
from .core import best_domains
from .io import atomic_open, read_json, read_tsv, sha256_file, write_json

log = logging.getLogger(__name__)

SUPPLEMENT = [
    ("report.md", "Database audit: superfamilies, tiers, entry groups at each identity, label curation, conflicts"),
    ("families.tsv", "Entries, entry groups, species, tiers per superfamily and family"),
    ("region_support.tsv", "Sequences and identity clusters per family and region"),
    ("label_conflicts.tsv", "Near-identical sequences with different superfamily labels"),
    ("unmapped_labels.tsv", "UniProt family labels without a curated mapping"),
    ("partition.md", "Test set, cross-validation folds, exclusions"),
    ("partition_summary.tsv", "Groups and entries per superfamily in train / test"),
    ("negatives.tsv", "Every negative sequence: source, role, QC removals and reasons"),
    ("profiles.tsv", "Every profile built: set, unit, sequences, match states"),
    ("cutoffs.tsv", "Baseline per-profile HMM cutoffs"),
    ("selection.md", "Choice of the defaults by grouped cross-validation on training rows (rule, intervals)"),
    ("selection_classification.tsv", "Cross-validated classifier variants at every tau"),
    ("selection_regions.tsv", "Cross-validated region-prediction variants"),
    ("selection_region_errors.tsv", "Boundary errors of region prediction in training cross-validation"),
    ("selection_recommended.json", "Settings recommended by cross-validation"),
    ("fit.md", "Classifier settings, posterior summary, held-out results"),
    ("fit/bayes_summary.csv", "Full posterior summary of the Bayesian model (R-hat, ESS)"),
    ("fit/bayes_status.json", "Whether the Bayesian posterior was sampled, reused or adopted"),
    ("predictions.tsv", "Held-out predictions of every classifier, with P(NONE) and best family"),
    ("evaluation.md", "Held-out performance by sequence type, identity, calibration, misses, ablations"),
    ("evaluation_metrics.tsv", "All evaluation metrics (long format)"),
    ("evaluation_calibration.tsv", "Calibration bins"),
    ("evaluation_misses.tsv", "Every misclassified held-out positive"),
    ("evaluation_ablations.tsv", "Test-set results of every classifier variant, with intervals"),
    ("anatomy.md", "Region prediction benchmark on the test set, every region variant"),
    ("anatomy_metrics.tsv", "Region prediction metrics per region, variant, family source, subset"),
    ("anatomy_sequences.tsv", "Per-sequence region prediction metrics (mature peptide, residue accuracy)"),
    ("anatomy_errors.tsv", "Residues between predicted and true boundaries, per variant"),
    ("anatomy_predictions.tsv", "Every predicted and true region boundary"),
    ("manifest.json", "Build manifest: versions, settings, UniProt release, file checksums"),
]


def _quantile_cut(values: list[float], q: float) -> float:
    return float(np.quantile(np.array(values), q)) if values else float("inf")


def calibrate_hmms(src: Path, dst: Path, negatives: list[tuple[str, str]], fpr: float, sc) -> list[dict]:
    from pyhmmer import plan7
    if not src.exists() or src.stat().st_size == 0:
        return []
    with plan7.HMMFile(src) as fh:
        hmms = list(fh)
    hits = best_domains(hmms, negatives, min_bits=sc.hmm_min_bits, max_mode=sc.max_mode, cpus=sc.threads)
    q = 1.0 - fpr / max(1, len(hmms))
    rows = []
    with atomic_open(dst, "wb") as fh:
        for h in hmms:
            name = h.name if isinstance(h.name, str) else h.name.decode()
            vals = [v.score for (sid, n), v in hits.items() if n == name]
            vals += [sc.hmm_min_bits] * (len(negatives) - len(vals))   # no domain reported
            cut = max(0.0, round(_quantile_cut(vals, q), 2))
            h.cutoffs.gathering = (cut, cut)
            h.write(fh)
            rows.append({"profile": name, "method": "hmm", "cutoff": cut, "negatives": len(negatives),
                         "quantile": round(q, 6)})
    return rows


def calibrate_pssms(src: Path, dst: Path, negatives: list[tuple[str, str]], fpr: float, pfsearch: str | None) -> list[dict]:
    """Rewrite /CUT_OFF (level 0) of each PROSITE profile from negatives' raw scores."""
    if not src.exists():
        return []
    text = src.read_text()
    blocks = [b + "//\n" for b in text.split("//\n") if b.strip()]
    if pfsearch is None:
        shutil.copyfile(src, dst)
        return [{"profile": "*", "method": "pssm", "cutoff": "uncalibrated (pfsearch not found)"}]
    from .score import pssm_scores
    q = 1.0 - fpr / max(1, len(blocks))
    seqs = [{"id": i, "sequence": s} for i, s in negatives]
    rows, out = [], []
    for blk in blocks:
        name = re.search(r"^ID   (\S+);", blk, re.M).group(1)
        tmp = dst.parent / f".{name}.prf"
        tmp.write_text(blk)
        try:
            vals = [r["score"] for r in pssm_scores(tmp, name, seqs, "package", pfsearch, [])]
        finally:
            tmp.unlink(missing_ok=True)
        cut = int(round(_quantile_cut(vals, q)))
        m = re.search(r"/NORMALIZATION:.*?R1=([-\d.]+);\s*R2=([-\d.]+);", blk)
        r1, r2 = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.01)
        lines = [ln for ln in blk.splitlines() if "/CUT_OFF:" not in ln]
        i = next((k for k, ln in enumerate(lines) if "/NORMALIZATION:" in ln), len(lines) - 1)
        lines.insert(i + 1, f"MA   /CUT_OFF: LEVEL=0; SCORE={cut}; N_SCORE={r1 + r2 * cut:.8f}; MODE=1; TEXT='!';")
        out.append("\n".join(lines) + "\n")
        rows.append({"profile": name, "method": "pssm", "cutoff": cut, "negatives": len(negatives), "quantile": round(q, 6)})
    with atomic_open(dst) as fh:
        fh.write("".join(out))
    return rows


def framework_tables(regions: list[dict], partition: list[dict]) -> dict:
    part = {p["accession"]: p for p in partition}
    counts: dict[str, Counter] = defaultdict(Counter)
    for r in regions:
        p = part.get(r["accession"])
        if r["region"] == "MAT" and p and p["scope"] == "conopeptide" and p["tier"] in ("1", "2") \
                and p["role"] in ("train", "test", "temporal_test"):
            counts[safe(p["superfamily"])][r["framework"] or "-"] += 1
    vocab = {fw for c in counts.values() for fw in c}
    return {"families": {f: dict(sorted(c.items())) for f, c in sorted(counts.items())}, "vocab_size": len(vocab) + 1}


def run_package(cfg: Config, rdir: Path, release: str, config_path: Path | None) -> dict:
    import shutil as _sh

    pc = cfg.package
    name = f"conodictor-db-{release}"
    pkg = rdir / "package" / name
    if pkg.exists():
        _sh.rmtree(pkg)
    pkg.mkdir(parents=True)
    prof = rdir / "profiles"

    negs = [(n["id"], n["sequence"]) for n in read_tsv(rdir / "negatives.tsv") if n["role"] in ("train", "test")]
    cut_rows = calibrate_hmms(prof / "ship.hmm", pkg / "profiles.hmm", negs, pc.cutoff_fpr, cfg.score)
    pfsearch = shutil.which(cfg.score.pfsearch)
    if (prof / "ship.prf").exists():
        cut_rows += calibrate_pssms(prof / "ship.prf", pkg / "profiles.prf", negs, pc.cutoff_fpr, pfsearch)
    for src, dst in (("ship.pwm.json", "profiles.pwm.json"), ("ship.full.hmm", "precursors.hmm"),
                     ("ship.full.labels.json", "precursors.labels.json"), ("ship.sites.json", "sites.json")):
        if (prof / src).exists():
            shutil.copyfile(prof / src, pkg / dst)
    write_json(pkg / "framework.json", framework_tables(read_tsv(rdir / "regions.tsv"), read_tsv(rdir / "partition.tsv")))
    model = read_json(rdir / "fit" / "model_ship.json") if (rdir / "fit" / "model_ship.json").exists() else None
    if not model or model.get("model") is None:
        model = read_json(rdir / "fit" / "model_map.json")
    write_json(pkg / "model.json", model)
    write_json(pkg / "settings.json", {s: cfg.section(s) for s in ("score", "fit", "anatomy", "build", "package")})
    with atomic_open(pkg / "cutoffs.tsv") as fh:
        fields = ["profile", "method", "cutoff", "negatives", "quantile"]
        fh.write("\t".join(fields) + "\n")
        for r in cut_rows:
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")

    # supplementary bundle
    sup = pkg / "supplementary"
    sup.mkdir()
    index = []                     # numbers are fixed per table (stable across releases); absent tables leave a gap
    for i, (rel, desc) in enumerate(SUPPLEMENT, 1):
        p = rdir / rel
        if p.exists():
            target = sup / f"S{i:02d}_{Path(rel).name}"
            shutil.copyfile(p, target)
            index.append((target.name, desc))
    num = len(SUPPLEMENT)
    for extra in ("labels", "entry_overrides"):
        p = getattr(cfg.labels, "overrides" if extra == "labels" else "entry_overrides")
        if p is not None and Path(p).exists():
            num += 1
            target = sup / f"S{num:02d}_{Path(p).name}"
            shutil.copyfile(p, target)
            index.append((target.name, "Curated label mapping" if extra == "labels" else "Per-entry label corrections"))
    if config_path is not None and Path(config_path).exists():
        target = sup / f"S{num + 1:02d}_{Path(config_path).name}"
        shutil.copyfile(config_path, target)
        index.append((target.name, "Configuration used for this build"))
    with atomic_open(sup / "README.md") as fh:
        fh.write(f"# Supplementary material: ConoDictor database, UniProt {release}\n\n")
        fh.write("| File | Content |\n|---|---|\n")
        fh.writelines(f"| `{f}` | {d} |\n" for f, d in index)

    import pyhmmer
    files = {p.relative_to(pkg).as_posix(): sha256_file(p) for p in sorted(pkg.rglob("*")) if p.is_file()}
    manifest = {
        "database": name, "uniprot_release": release, "conodb_version": __version__,
        "python": sys.version.split()[0], "platform": platform.platform(),
        "packages": {"pyhmmer": pyhmmer.__version__}, "families": model.get("families"),
        "model": model.get("model"), "features": model.get("features"),
        "pssm_cutoffs": "calibrated" if pfsearch else "uncalibrated", "files": files,
    }
    write_json(pkg / "manifest.json", manifest)
    with atomic_open(pkg / "README.md") as fh:
        fh.write(README.format(name=name, release=release, version=__version__,
                               families=", ".join(model.get("families") or []),
                               features=", ".join(model.get("features") or [])))
    archive = None
    if pc.archive:
        archive = rdir / "package" / f"{name}.tar.gz"
        with tarfile.open(archive, "w:gz") as tf:
            tf.add(pkg, arcname=name)
    return {"package": str(pkg.relative_to(rdir)), "files": len(files) + 2,
            "archive": archive and str(archive.relative_to(rdir)), "profiles_calibrated": len(cut_rows)}


README = """# {name}

ConoDictor reference database built by conodb {version} from UniProt release {release}.

| File | Content |
|---|---|
| `profiles.hmm` | Region HMMs (SIG/PRO/MAT/POST) per superfamily; GA cutoffs calibrated on negatives |
| `profiles.prf` | Region PSSMs (pftools generalized profiles); /CUT_OFF calibrated on negatives |
| `profiles.pwm.json` | Ungapped lightmotif windows per region |
| `precursors.hmm` + `precursors.labels.json` | Full-precursor HMMs and the region of each match state |
| `sites.json` | Cleavage-site windows (signal peptide end, mature start, mature end) |
| `framework.json` | Cysteine-framework counts per family |
| `model.json` | Classifier: conditional logit over {families} + NONE; features {features} |
| `settings.json` | Settings needed to compute features exactly as in training |
| `cutoffs.tsv` | The calibrated cutoffs |
| `manifest.json` | Versions and sha256 of every file |
| `supplementary/` | Tables behind the evaluation (index in its README) |

Prediction must use conodb.core (shared with conodictor) so that features are
computed exactly as during training.
"""
