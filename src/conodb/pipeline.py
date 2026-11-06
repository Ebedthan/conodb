"""Stage runner: one command rebuilds the database, redoing only what changed.

Layout of the output root:

    db/
      LATEST                      name of the last successfully built release
      latest -> releases/<rel>    convenience symlink (where supported)
      releases/<rel>/
        raw/uniprot.jsonl.gz      exact UniProt download
        entries.tsv  regions.tsv  regions/<REGION>.fasta
        unmapped_labels.tsv  parse_warnings.tsv
        clusters.tsv  entry_groups.tsv
        families.tsv  region_support.tsv  label_conflicts.tsv
        diff_entries.tsv  diff_tiers.tsv  report.md
        manifest.json             versions, config, hashes of every output
        .stamps/<stage>.json      input hash used to skip unchanged stages
        logs/<timestamp>.log

A stage is skipped when its key (hash of: stage name, conodb version, the
config sections it reads, and the hashes of its input files) matches the
stamp from its last successful run and its outputs are intact. Adding a
stage = writing a function and appending a Stage to STAGES.
"""

from __future__ import annotations

import logging
import os
import platform
import shutil
import sys
import time
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .audit import (
    CONFLICT_FIELDS, DIFF_FIELDS, diff_entries, family_fields, diff_tiers, render_report,
    run_audit, support_fields,
)
from .cluster import cluster_levels
from .config import Config
from .io import (
    atomic_open, iter_jsonl_gz, read_json, read_tsv, sha256_file, sha256_obj, write_fasta,
    write_json, write_tsv,
)
from .build import PROFILE_FIELDS, run_build
from .negatives import NEGATIVE_FIELDS, build_negatives, fetch_negatives
from .evaluate import ABLATION_FIELDS, CALIB_FIELDS, METRIC_FIELDS, MISS_FIELDS, render_evaluation, run_evaluate
from .anatomy import METRIC_FIELDS as ANATOMY_METRIC_FIELDS
from .anatomy import PRED_FIELDS as ANATOMY_PRED_FIELDS
from .anatomy import SEQ_FIELDS as ANATOMY_SEQ_FIELDS
from .anatomy import ERROR_FIELDS as ANATOMY_ERROR_FIELDS
from .anatomy import render_anatomy, run_anatomy
from .select import CLS_FIELDS as SELECT_CLS_FIELDS
from .select import ERR_FIELDS as SELECT_ERR_FIELDS
from .select import REG_FIELDS as SELECT_REG_FIELDS
from .select import render_selection, run_select
from .fit import (
    FIT_FORMAT, STAN_FILE, NothingToFit, bayes_core_key, bayes_wanted, load_bayes, load_design, render_fit, run_bayes,
    run_fit,
)
from .package import run_package
from .score import SCORE_FIELDS, run_score
from .partition import PARTITION_FIELDS, SUMMARY_FIELDS, render_partition, run_partition
from .parse import (
    ENTRY_FIELDS, PROFILE_REGIONS, REGION_FIELDS, LabelResolver, parse_entry, unmapped_rows,
)
from .uniprot import ReleaseInfo, fetch_all, probe

log = logging.getLogger("conodb")

EVAL_FORMAT = 3
ANATOMY_FORMAT = 3
SELECT_FORMAT = 2
PACKAGE_FORMAT = 1
ANATOMY_OUTPUTS = ("anatomy.md", "anatomy_metrics.tsv", "anatomy_sequences.tsv", "anatomy_predictions.tsv",
                   "anatomy_errors.tsv")
SELECT_OUTPUTS = ("selection.md", "selection_classification.tsv", "selection_regions.tsv", "selection_region_errors.tsv",
                  "selection_recommended.json")


class CurationNeeded(RuntimeError):
    """Raised in strict mode when UniProt has family labels nobody curated."""


@dataclass
class Context:
    cfg: Config
    root: Path
    release: ReleaseInfo
    rdir: Path
    local_input: Path | None = None
    previous: Path | None = None
    stats: dict = field(default_factory=dict)

    def p(self, rel: str) -> Path:
        return self.rdir / rel


@dataclass
class Stage:
    name: str
    sections: tuple[str, ...]
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    func: Callable[[Context], dict]
    extra_key: Callable[[Context], object] = lambda ctx: None
    verify: Callable[[Context], bool] | None = None
    description: str = ""
    # Replaces the generic key (used by the expensive Bayesian stage, whose
    # key must not change with the conodb version or unrelated settings).
    key_fn: Callable[[Context], str] | None = None


# ------------------------------------------------------------------ stages --

def stage_fetch(ctx: Context) -> dict:
    out = ctx.p("raw/uniprot.jsonl.gz")
    if ctx.local_input is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ctx.local_input, out)
        n = sum(1 for _ in iter_jsonl_gz(out))
        ctx.release.total = n
        return {"entries": n, "source": str(ctx.local_input)}
    info = fetch_all(ctx.cfg.source, out)
    if info.release != ctx.release.release:
        raise RuntimeError(f"UniProt release changed from {ctx.release.release} to {info.release}; rerun")
    ctx.release.total = info.total
    return {"entries": info.total}


def stage_parse(ctx: Context) -> dict:
    resolver = LabelResolver.from_config(ctx.cfg.labels)
    entries, regions, warnings = [], [], []
    for raw in iter_jsonl_gz(ctx.p("raw/uniprot.jsonl.gz")):
        e, r, w = parse_entry(raw, resolver, ctx.cfg.regions)
        entries.append(e)
        regions.extend(r)
        warnings.extend(w)
    entries.sort(key=lambda e: e["accession"])
    regions.sort(key=lambda r: (r["accession"], r["start"]))
    write_tsv(ctx.p("entries.tsv"), entries, ENTRY_FIELDS)
    write_tsv(ctx.p("regions.tsv"), regions, REGION_FIELDS)
    unmapped = unmapped_rows(resolver)
    write_tsv(ctx.p("unmapped_labels.tsv"), unmapped, ["raw_superfamily", "raw_family", "entries"])
    write_tsv(ctx.p("parse_warnings.tsv"), ({"warning": w} for w in warnings), ["warning"])
    for region in PROFILE_REGIONS:
        write_fasta(
            ctx.p(f"regions/{region}.fasta"),
            ((r["region_id"], r["sequence"]) for r in regions if r["region"] == region),
        )
    if unmapped and ctx.cfg.labels.strict:
        raise CurationNeeded(
            f"{len(unmapped)} UniProt family label(s) need curation; see {ctx.p('unmapped_labels.tsv')} "
            f"and add them to {ctx.cfg.labels.overrides or 'a labels.tsv file'}"
        )
    scopes = Counter(e["scope"] for e in entries)
    return {"entries": len(entries), "regions": len(regions), "warnings": len(warnings),
            "unmapped_labels": len(unmapped), "scopes": dict(scopes)}


def stage_cluster(ctx: Context) -> dict:
    entries = read_tsv(ctx.p("entries.tsv"))
    regions = read_tsv(ctx.p("regions.tsv"))
    idents = ctx.cfg.cluster.identities
    threads = ctx.cfg.cluster.threads
    rows = []
    sets = {"FULL": [(e["accession"], e["sequence"]) for e in entries if e["scope"] != "exclude" and e["sequence"]]}
    keep = {e["accession"] for e in entries if e["scope"] != "exclude"}
    for region in PROFILE_REGIONS:
        sets[region] = [(r["region_id"], r["sequence"]) for r in regions if r["region"] == region and r["accession"] in keep]
    counts = {}
    for kind, items in sets.items():
        t0 = time.time()
        levels = cluster_levels(items, idents, threads)
        for t, assign in levels.items():
            rows.extend({"item_id": i, "kind": kind, "identity": t, "cluster": c} for i, c in assign.items())
        counts[kind] = {f"{t:g}": len(set(a.values())) for t, a in levels.items()}
        log.info("  %-4s %5d sequences -> %s clusters (%.1fs)", kind, len(items),
                 "/".join(str(counts[kind][f"{t:g}"]) for t in sorted(idents, reverse=True)), time.time() - t0)
    rows.sort(key=lambda r: (r["kind"], -r["identity"], r["item_id"]))
    write_tsv(ctx.p("clusters.tsv"), rows, ["item_id", "kind", "identity", "cluster"])
    g = ctx.cfg.cluster.group_identity
    groups = [{"accession": r["item_id"], "group": r["cluster"]} for r in rows if r["kind"] == "FULL" and r["identity"] == g]
    write_tsv(ctx.p("entry_groups.tsv"), groups, ["accession", "group"])
    return {"clusters": counts}


def stage_audit(ctx: Context) -> dict:
    cfg = ctx.cfg
    entries = read_tsv(ctx.p("entries.tsv"))
    regions = read_tsv(ctx.p("regions.tsv"))
    clusters = read_tsv(ctx.p("clusters.tsv"))
    groups = {r["accession"]: r["group"] for r in read_tsv(ctx.p("entry_groups.tsv"))}
    result = run_audit(cfg, entries, regions, clusters, groups)
    write_tsv(ctx.p("families.tsv"), result.families, family_fields(cfg))
    write_tsv(ctx.p("region_support.tsv"), result.region_support, support_fields(cfg))
    write_tsv(ctx.p("label_conflicts.tsv"), result.conflicts, CONFLICT_FIELDS)

    warn_counts = Counter()
    for w in read_tsv(ctx.p("parse_warnings.tsv")):
        msg = w["warning"].split(": ", 1)[-1]
        warn_counts[_warning_kind(msg)] += 1

    entry_diff = tier_diff = None
    prev_name = None
    if ctx.previous is not None:
        prev_name = ctx.previous.name
        entry_diff = diff_entries(read_tsv(ctx.previous / "entries.tsv"), entries)
        tier_diff = diff_tiers(read_tsv(ctx.previous / "families.tsv"), [
            {k: str(v) for k, v in f.items()} for f in result.families])
    write_tsv(ctx.p("diff_entries.tsv"), entry_diff or [], DIFF_FIELDS)
    write_tsv(ctx.p("diff_tiers.tsv"), tier_diff or [], ["superfamily", "family", "old", "new"])

    report = render_report(
        cfg, {"release": ctx.release.release, "release_date": ctx.release.release_date, "total": len(entries)},
        result, read_tsv(ctx.p("unmapped_labels.tsv")), dict(warn_counts), entry_diff, tier_diff, prev_name,
    )
    with atomic_open(ctx.p("report.md")) as fh:
        fh.write(report)
    tiers = Counter(f["tier"] for f in result.families if f["family"] == "*" and f["tier"] != "")
    return {"tiers": {str(k): v for k, v in sorted(tiers.items())}, "label_conflicts": len(result.conflicts),
            "changes": len(entry_diff) if entry_diff is not None else None}


def stage_partition(ctx: Context) -> dict:
    entries = read_tsv(ctx.p("entries.tsv"))
    families = read_tsv(ctx.p("families.tsv"))
    audit_tiers = {f["superfamily"]: f["tier"] for f in families if f["family"] == "*" and f["scope"] == "conopeptide"}
    prev_file = ctx.previous / "partition.tsv" if ctx.previous else None
    previous = read_tsv(prev_file) if prev_file and prev_file.exists() else None
    res = run_partition(ctx.cfg, entries, read_tsv(ctx.p("entry_groups.tsv")), read_tsv(ctx.p("label_conflicts.tsv")),
                        audit_tiers, previous)
    write_tsv(ctx.p("partition.tsv"), res.rows, PARTITION_FIELDS)
    write_tsv(ctx.p("partition_summary.tsv"), res.summary, SUMMARY_FIELDS)
    with atomic_open(ctx.p("partition.md")) as fh:
        fh.write(render_partition(ctx.cfg, ctx.release.release, res, ctx.previous.name if previous else None))
    roles = Counter(r["role"] for r in res.rows)
    return {"roles": dict(roles), "excluded": len(res.excluded), "moved_to_test": len(res.reassigned)}


def stage_negatives(ctx: Context) -> dict:
    nc = ctx.cfg.negatives
    raw = ctx.p("raw/negatives.jsonl.gz")
    fetched = {}
    if nc.enabled:
        fetched = fetch_negatives(ctx.cfg, ctx.release, raw, offline=ctx.local_input is not None)
    prev_file = ctx.previous / "negatives.tsv" if ctx.previous else None
    previous = read_tsv(prev_file) if prev_file and prev_file.exists() else None
    rows, stats = build_negatives(ctx.cfg, read_tsv(ctx.p("entries.tsv")), read_tsv(ctx.p("regions.tsv")),
                                  read_tsv(ctx.p("partition.tsv")), raw if nc.enabled else None, previous)
    write_tsv(ctx.p("negatives.tsv"), rows, NEGATIVE_FIELDS)
    write_fasta(ctx.p("negatives.fasta"), ((r["id"], r["sequence"]) for r in rows if r["role"] in ("train", "test")))
    return stats | {"fetch": fetched}


def stage_build(ctx: Context) -> dict:
    rows, stats = run_build(ctx.cfg, ctx.rdir, read_tsv(ctx.p("entries.tsv")), read_tsv(ctx.p("regions.tsv")),
                            read_tsv(ctx.p("partition.tsv")))
    write_tsv(ctx.p("profiles.tsv"), rows, PROFILE_FIELDS)
    for note in stats.get("notes", []):
        log.info("  %s", note)
    return stats


def stage_score(ctx: Context) -> dict:
    scores, features, fields, stats = run_score(
        ctx.cfg, ctx.rdir, read_tsv(ctx.p("entries.tsv")), read_tsv(ctx.p("partition.tsv")),
        read_tsv(ctx.p("negatives.tsv")), read_tsv(ctx.p("profiles.tsv")), read_tsv(ctx.p("regions.tsv")))
    write_tsv(ctx.p("scores.tsv"), scores, SCORE_FIELDS)
    write_tsv(ctx.p("features.tsv"), features, fields)
    return stats


# The former single 'fit' stage (conodb 0.1.0, fit format 4) sampled the
# Bayesian model itself. Its key is recomputed here so a posterior it
# produced is adopted instead of being sampled again.
_LEGACY_FIT_FORMAT = 4


def _legacy_fit_key(ctx: Context) -> str:
    return sha256_obj({
        "stage": "fit",
        "conodb": "0.1.0",
        "config": {s: ctx.cfg.section(s) for s in ("fit", "score")},
        "inputs": {"features.tsv": sha256_file(ctx.p("features.tsv"))},
        "extra": [ctx.cfg.fit.bayes != "false" and _cmdstan_id(ctx.cfg), _LEGACY_FIT_FORMAT, sha256_file(STAN_FILE)],
    })


def _legacy_match(ctx: Context) -> bool:
    stamp = ctx.p(".stamps/fit.json")
    if not stamp.exists():
        return False
    try:
        return read_json(stamp).get("key") == _legacy_fit_key(ctx)
    except Exception:
        return False


def stage_bayes(ctx: Context) -> dict:
    try:
        status = run_bayes(ctx.cfg, ctx.rdir, ctx.root / ".cache" / "stan", _legacy_match(ctx))
    except NothingToFit as e:
        status = {"status": "skipped", "reason": str(e)}
    write_json(ctx.p("fit/bayes_status.json"), status)
    if status["status"] == "skipped":
        log.info("  Bayesian fit skipped: %s", status.get("reason"))
    return {k: v for k, v in status.items() if k != "key"}


def _bayes_key(ctx: Context) -> str:
    try:
        core = bayes_core_key(ctx.cfg, load_design(ctx.cfg, ctx.rdir))
    except NothingToFit:
        core = None
    return sha256_obj({"core": core, "wanted": bayes_wanted(ctx.cfg)})


def _verify_bayes(ctx: Context) -> bool:
    status = read_json(ctx.p("fit/bayes_status.json")).get("status")
    return status == "skipped" or load_bayes(ctx.cfg, ctx.rdir) is not None


def stage_fit(ctx: Context) -> dict:
    try:
        out = run_fit(ctx.cfg, ctx.rdir)
    except NothingToFit as e:
        log.warning("Nothing to fit: %s", e)
        write_tsv(ctx.p("cutoffs.tsv"), [], ["family", "region", "method", "cutoff_bits", "negatives", "quantile"])
        write_tsv(ctx.p("predictions.tsv"), [], ["id"])
        write_json(ctx.p("fit/model_map.json"), {"model": None, "reason": str(e)})
        write_json(ctx.p("fit/model_ship.json"), {"model": None, "reason": str(e)})
        with atomic_open(ctx.p("fit.md")) as fh:
            fh.write(f"# Classifier fit: UniProt {ctx.release.release}\n\nNothing to fit: {e}.\n")
        return {"nothing_to_fit": str(e)}
    write_tsv(ctx.p("cutoffs.tsv"), out["cutoffs"], ["family", "region", "method", "cutoff_bits", "negatives", "quantile"])
    write_tsv(ctx.p("predictions.tsv"), out["predictions"], out["pred_fields"])
    if not ctx.cfg.fit.ship_all_data:
        write_json(ctx.p("fit/model_ship.json"), {"model": None, "reason": "fit.ship_all_data = false"})
    with atomic_open(ctx.p("fit.md")) as fh:
        fh.write(render_fit(ctx.cfg, ctx.release.release, out))
    summary = {split: {m: {k: round(v, 4) for k, v in r.items() if isinstance(v, float)}
                       for m, r in res.items() if r} for split, res in out["results"].items()}
    return {"results": summary}


def stage_anatomy(ctx: Context) -> dict:
    if not ctx.cfg.build.full_models:
        for o in ANATOMY_OUTPUTS[1:]:
            write_tsv(ctx.p(o), [], ["empty"])
        with atomic_open(ctx.p("anatomy.md")) as fh:
            fh.write("# Region prediction benchmark\n\nbuild.full_models = false.\n")
        return {"evaluated": False}
    out = run_anatomy(ctx.cfg, ctx.rdir)
    write_tsv(ctx.p("anatomy_predictions.tsv"), out["predictions"], ANATOMY_PRED_FIELDS)
    write_tsv(ctx.p("anatomy_metrics.tsv"), out["metrics"], ANATOMY_METRIC_FIELDS)
    write_tsv(ctx.p("anatomy_sequences.tsv"), out["sequences"], ANATOMY_SEQ_FIELDS)
    write_tsv(ctx.p("anatomy_errors.tsv"), out["errors"], ANATOMY_ERROR_FIELDS)
    rec = read_json(ctx.p("selection_recommended.json")) if ctx.p("selection_recommended.json").exists() else {}
    with atomic_open(ctx.p("anatomy.md")) as fh:
        fh.write(render_anatomy(ctx.cfg, ctx.release.release, out, rec.get("region_variant")))
    return {"sequences": out["n"], "predictions": len(out["predictions"])}


def stage_select(ctx: Context) -> dict:
    if ctx.cfg.evaluate.select:
        try:
            out = run_select(ctx.cfg, ctx.rdir)
        except NothingToFit as e:
            out, reason = None, str(e)
    else:
        out, reason = None, "evaluate.select = false"
    if out is None:
        write_tsv(ctx.p("selection_classification.tsv"), [], ["empty"])
        write_tsv(ctx.p("selection_regions.tsv"), [], ["empty"])
        write_tsv(ctx.p("selection_region_errors.tsv"), [], ["empty"])
        write_json(ctx.p("selection_recommended.json"), {})
        with atomic_open(ctx.p("selection.md")) as fh:
            fh.write(f"# Selection of the defaults\n\nNot run: {reason}.\n")
        return {"selected": False, "reason": reason}
    write_tsv(ctx.p("selection_classification.tsv"), out["classifier"]["table"], SELECT_CLS_FIELDS)
    write_tsv(ctx.p("selection_regions.tsv"), out["regions"]["table"], SELECT_REG_FIELDS)
    write_tsv(ctx.p("selection_region_errors.tsv"), out["regions"].get("errors", []), SELECT_ERR_FIELDS)
    write_json(ctx.p("selection_recommended.json"), out["recommended"] | {"changes": out["changes"]})
    with atomic_open(ctx.p("selection.md")) as fh:
        fh.write(render_selection(ctx.cfg, ctx.release.release, out))
    if out["changes"]:
        log.info("  cross-validation recommends: %s (see selection.md)", "; ".join(out["changes"]))
    return {"classifier": out["recommended"]["classifier_variant"],
            "regions": out["recommended"].get("region_variant"), "changes": out["changes"]}


def stage_package(ctx: Context) -> dict:
    if read_json(ctx.p("fit/model_map.json")).get("model") is None:
        write_json(ctx.p("package.json"), {"package": ".", "skipped": "no model"})
        return {"packaged": False}
    out = run_package(ctx.cfg, ctx.rdir, ctx.release.release, ctx.cfg.path)
    write_json(ctx.p("package.json"), out)
    return out


EVAL_OUTPUTS = ("evaluation.md", "evaluation_metrics.tsv", "evaluation_misses.tsv", "evaluation_calibration.tsv",
                "evaluation_ablations.tsv")


def stage_evaluate(ctx: Context) -> dict:
    if read_json(ctx.p("fit/model_map.json")).get("model") is None:
        for o in EVAL_OUTPUTS[1:]:
            write_tsv(ctx.p(o), [], ["empty"])
        with atomic_open(ctx.p("evaluation.md")) as fh:
            fh.write(f"# Evaluation: UniProt {ctx.release.release}\n\nNo model to evaluate.\n")
        return {"evaluated": False}
    out = run_evaluate(ctx.cfg, ctx.rdir)
    write_tsv(ctx.p("evaluation_metrics.tsv"), out["metrics"], METRIC_FIELDS)
    write_tsv(ctx.p("evaluation_misses.tsv"), out["misses"], MISS_FIELDS)
    write_tsv(ctx.p("evaluation_calibration.tsv"), out["calibration"], CALIB_FIELDS)
    write_tsv(ctx.p("evaluation_ablations.tsv"), out["ablations"], ABLATION_FIELDS)
    with atomic_open(ctx.p("evaluation.md")) as fh:
        fh.write(render_evaluation(ctx.cfg, ctx.release.release, out))
    return {"methods": out["methods"], "misses": len(out["misses"]), "ablations": len(out["ablations"])}


def _verify_profiles(ctx: Context) -> bool:
    """Every profile file listed in profiles.tsv still exists unchanged."""
    for r in read_tsv(ctx.p("profiles.tsv")):
        for col, sha in (("hmm", "hmm_sha256"), ("pssm", "pssm_sha256")):
            if r.get(col) and (not ctx.p(r[col]).exists() or sha256_file(ctx.p(r[col])) != r[sha]):
                return False
    return True


def _warning_kind(msg: str) -> str:
    for key in ("uncertain boundary", "unknown boundary", "outside sequence", "shorter than"):
        if key in msg:
            return msg.split(" ", 1)[0] + " " + key
    return msg


STAGES: list[Stage] = [
    Stage("fetch", ("source",), (), ("raw/uniprot.jsonl.gz",), stage_fetch,
          extra_key=lambda ctx: [ctx.release.release, str(ctx.local_input and sha256_file(ctx.local_input))],
          description="Download matching UniProt entries (JSON)"),
    Stage("parse", ("labels", "regions"), ("raw/uniprot.jsonl.gz",),
          ("entries.tsv", "regions.tsv", "unmapped_labels.tsv", "parse_warnings.tsv")
          + tuple(f"regions/{r}.fasta" for r in PROFILE_REGIONS),
          stage_parse,
          extra_key=lambda ctx: [p and sha256_file(p) for p in (ctx.cfg.labels.overrides, ctx.cfg.labels.entry_overrides)],
          description="Labels, precursor regions, evidence"),
    Stage("cluster", ("cluster",), ("entries.tsv", "regions.tsv"), ("clusters.tsv", "entry_groups.tsv"),
          stage_cluster, description="Identity clusters per region and entry groups"),
    Stage("audit", ("tiers", "cluster", "audit"),
          ("entries.tsv", "regions.tsv", "clusters.tsv", "entry_groups.tsv", "unmapped_labels.tsv", "parse_warnings.tsv"),
          ("families.tsv", "region_support.tsv", "label_conflicts.tsv", "diff_entries.tsv", "diff_tiers.tsv", "report.md"),
          stage_audit,
          extra_key=lambda ctx: ctx.previous and sha256_file(ctx.previous / "entries.tsv"),
          description="Support per family/region, tiers, QC, diff, report"),
    Stage("partition", ("partition", "tiers", "cluster"),
          ("entries.tsv", "entry_groups.tsv", "families.tsv", "label_conflicts.tsv"),
          ("partition.tsv", "partition_summary.tsv", "partition.md"),
          stage_partition,
          extra_key=lambda ctx: ctx.previous and (ctx.previous / "partition.tsv").exists()
          and sha256_file(ctx.previous / "partition.tsv"),
          description="Held-out test set and CV folds (stable across releases)"),
    Stage("negatives", ("negatives", "partition", "cluster"), ("entries.tsv", "regions.tsv", "partition.tsv"),
          ("negatives.tsv", "negatives.fasta"), stage_negatives,
          extra_key=lambda ctx: [ctx.release.release, ctx.local_input is not None,
                                 ctx.previous and (ctx.previous / "negatives.tsv").exists()
                                 and sha256_file(ctx.previous / "negatives.tsv")],
          description="Negatives: UniProt query, venom 'other', shuffled precursors"),
    Stage("build", ("build", "partition", "anatomy"), ("entries.tsv", "regions.tsv", "partition.tsv"),
          ("profiles.tsv", "profiles/final.hmm", "profiles/final.sites.json", "profiles/ship.sites.json"), stage_build,
          extra_key=lambda ctx: [shutil.which(ctx.cfg.build.pfw), shutil.which(ctx.cfg.build.pfmake)]
          if ctx.cfg.build.pssm != "false" else None,
          verify=_verify_profiles,
          description="Region HMM/PSSM/PWM profiles, full-precursor models, site windows: final, holdouts, ship"),
    Stage("score", ("score", "anatomy"), ("entries.tsv", "regions.tsv", "partition.tsv", "negatives.tsv", "profiles.tsv"),
          ("scores.tsv", "features.tsv"), stage_score,
          extra_key=lambda ctx: shutil.which(ctx.cfg.score.pfsearch) if ctx.cfg.score.pssm != "false" else None,
          description="Out-of-fold HMM/PSSM/PWM scores and framework feature of positives and negatives"),
    Stage("select", ("evaluate", "fit", "anatomy", "score", "partition"),
          ("entries.tsv", "regions.tsv", "partition.tsv", "features.tsv", "profiles.tsv"),
          SELECT_OUTPUTS, stage_select, extra_key=lambda ctx: [SELECT_FORMAT],
          description="Choice of the defaults: grouped cross-validation of classifier and region variants on training rows"),
    Stage("bayes", ("fit",), ("features.tsv",), ("fit/bayes_status.json",), stage_bayes,
          key_fn=_bayes_key, verify=_verify_bayes,
          description="Hierarchical Bayesian model (Stan); reused whenever the training design and model settings are unchanged"),
    Stage("fit", ("fit", "score"), ("features.tsv", "fit/bayes_status.json"),
          ("cutoffs.tsv", "predictions.tsv", "fit.md", "fit/model_map.json", "fit/model_ship.json"), stage_fit,
          extra_key=lambda ctx: [FIT_FORMAT],
          description="Baseline and penalised classifier (seconds); shipped model refitted on all data"),
    Stage("anatomy", ("anatomy", "score"), ("entries.tsv", "regions.tsv", "partition.tsv", "predictions.tsv",
                                            "profiles.tsv", "selection_recommended.json"),
          ANATOMY_OUTPUTS, stage_anatomy, extra_key=lambda ctx: [ANATOMY_FORMAT],
          description="Region prediction benchmark on the test set: every region variant"),
    Stage("evaluate", ("evaluate", "fit", "score", "cluster"),
          ("entries.tsv", "regions.tsv", "partition.tsv", "negatives.tsv", "features.tsv", "scores.tsv",
           "fit/model_map.json", "fit/bayes_status.json", "selection_recommended.json"),
          EVAL_OUTPUTS, stage_evaluate, extra_key=lambda ctx: [EVAL_FORMAT],
          description="Held-out performance by sequence type and identity, calibration, misses, ablations"),
    Stage("package", ("package", "score", "fit", "anatomy", "build"),
          ("profiles.tsv", "negatives.tsv", "regions.tsv", "partition.tsv", "fit/model_map.json",
           "fit/model_ship.json", "evaluation.md", "anatomy.md", "selection.md"),
          ("package.json",), stage_package, extra_key=lambda ctx: [PACKAGE_FORMAT, ctx.release.release],
          verify=lambda ctx: (ctx.rdir / read_json(ctx.p("package.json"))["package"]).exists(),
          description="Database for conodictor (calibrated cutoffs, models) + supplementary bundle"),
]


def _cmdstan_id(cfg) -> str | None:
    from .fit import cmdstan_available
    if not cmdstan_available(cfg):
        return None
    import cmdstanpy
    return cmdstanpy.cmdstan_path()
STAGE_NAMES = [s.name for s in STAGES]


# ------------------------------------------------------------------ runner --

def _stage_key(ctx: Context, stage: Stage) -> str:
    if stage.key_fn is not None:
        return stage.key_fn(ctx)
    return sha256_obj({
        "stage": stage.name,
        "conodb": __version__,
        "config": {s: ctx.cfg.section(s) for s in stage.sections},
        "inputs": {i: sha256_file(ctx.p(i)) for i in stage.inputs},
        "extra": stage.extra_key(ctx),
    })


def _up_to_date(ctx: Context, stage: Stage, key: str) -> bool:
    stamp = ctx.p(f".stamps/{stage.name}.json")
    if not stamp.exists():
        return False
    s = read_json(stamp)
    if s.get("key") != key:
        return False
    if not all(ctx.p(o).exists() and sha256_file(ctx.p(o)) == s["outputs"].get(o) for o in stage.outputs):
        return False
    return stage.verify(ctx) if stage.verify else True


def _find_previous(root: Path, current: Path) -> Path | None:
    best = None
    for d in (root / "releases").glob("*"):
        m = d / "manifest.json"
        if d == current or not m.exists():
            continue
        finished = read_json(m).get("finished", "")
        if best is None or finished > best[0]:
            best = (finished, d)
    return best[1] if best else None


def _setup_file_log(rdir: Path) -> logging.Handler:
    (rdir / "logs").mkdir(parents=True, exist_ok=True)
    h = logging.FileHandler(rdir / "logs" / f"{datetime.now().strftime('%Y%m%dT%H%M%S')}.log", encoding="utf-8")
    h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    h.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(h)
    return h


def resolve_release(cfg: Config, local_input: Path | None, release_name: str | None) -> ReleaseInfo:
    if local_input is not None:
        name = release_name or f"local-{sha256_file(local_input)[:10]}"
        return ReleaseInfo(release=name, release_date="", total=-1)
    return probe(cfg.source)


def run(cfg: Config, *, force: bool = False, from_stage: str | None = None, until_stage: str | None = None,
        local_input: Path | None = None, release_name: str | None = None) -> Path:
    root = cfg.output.root
    release = resolve_release(cfg, local_input, release_name)
    rdir = root / "releases" / release.dirname
    rdir.mkdir(parents=True, exist_ok=True)
    handler = _setup_file_log(rdir)
    ctx = Context(cfg=cfg, root=root, release=release, rdir=rdir, local_input=local_input,
                  previous=_find_previous(root, rdir))
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log.info("conodb %s: building release %s in %s", __version__, release.release, rdir)
    if ctx.previous:
        log.info("Previous release for comparison: %s", ctx.previous.name)

    names = STAGE_NAMES
    start = names.index(from_stage) if from_stage else 0
    stop = names.index(until_stage) + 1 if until_stage else len(names)
    ran, stage_info = [], {}
    try:
        for i, stage in enumerate(STAGES[:stop]):
            key = _stage_key(ctx, stage)
            forced = force or (from_stage is not None and i >= start)
            if not forced and _up_to_date(ctx, stage, key):
                log.info("[%s] up to date, skipped", stage.name)
                stage_info[stage.name] = read_json(ctx.p(f".stamps/{stage.name}.json"))
                continue
            log.info("[%s] %s", stage.name, stage.description)
            t0 = time.time()
            stats = stage.func(ctx)
            stamp = {
                "key": key,
                "outputs": {o: sha256_file(ctx.p(o)) for o in stage.outputs},
                "stats": stats,
                "seconds": round(time.time() - t0, 2),
                "finished": datetime.now(timezone.utc).isoformat(timespec="microseconds"),
            }
            write_json(ctx.p(f".stamps/{stage.name}.json"), stamp)
            stage_info[stage.name] = stamp
            ran.append(stage.name)
            log.info("[%s] done in %.1fs: %s", stage.name, stamp["seconds"], stats)
    finally:
        logging.getLogger().removeHandler(handler)
        handler.close()

    if stop == len(names):
        _write_manifest(ctx, started, stage_info)
        _update_latest(root, rdir)
    if not ran:
        log.info("Nothing to do: release %s is already built with this configuration.", release.release)
    else:
        log.info("Built release %s (stages run: %s). Report: %s", release.release, ", ".join(ran), ctx.p("report.md"))
    return rdir


def _write_manifest(ctx: Context, started: str, stage_info: dict) -> None:
    import pyfamsa
    import pyhmmer
    import pyopal

    outputs = sorted({o for s in STAGES for o in s.outputs})
    write_json(ctx.p("manifest.json"), {
        "conodb_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": {"pyopal": pyopal.__version__, "pyhmmer": pyhmmer.__version__, "pyfamsa": pyfamsa.__version__},
        "uniprot": {"release": ctx.release.release, "release_date": ctx.release.release_date,
                    "query": ctx.cfg.source.query,
                    "entries": stage_info.get("fetch", {}).get("stats", {}).get("entries", ctx.release.total)},
        "config": {s: ctx.cfg.section(s) for s in ("source", "labels", "regions", "cluster", "tiers", "audit", "partition", "negatives", "build", "score", "fit", "evaluate", "anatomy", "package")},
        "entry_overrides_sha256": ctx.cfg.labels.entry_overrides and sha256_file(ctx.cfg.labels.entry_overrides),
        "labels_sha256": ctx.cfg.labels.overrides and sha256_file(ctx.cfg.labels.overrides),
        "previous_release": ctx.previous.name if ctx.previous else None,
        "started": started,
        "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "stages": {k: {"stats": v.get("stats"), "finished": v.get("finished")} for k, v in stage_info.items()},
        "files": {o: sha256_file(ctx.p(o)) for o in outputs if ctx.p(o).exists()},
    })


GITIGNORE = """# Written by conodb. Commit the small summary tables of each release; the
# rest is regenerated by `conodb run` (or kept as a CI artifact / GitHub
# release asset for exact reproducibility).
latest
releases/*/raw/
releases/*/negatives.fasta
releases/*/scores.tsv
releases/*/fit/stan/
releases/*/package/
releases/*/profiles/ship/
releases/*/fit/bayes_draws.npz
releases/*/regions/
releases/*/clusters.tsv
releases/*/logs/
releases/*/.stamps/
releases/*/profiles/final/
releases/*/profiles/holdout/
"""


def _update_latest(root: Path, rdir: Path) -> None:
    if not (root / ".gitignore").exists():
        (root / ".gitignore").write_text(GITIGNORE, encoding="utf-8")
    with atomic_open(root / "LATEST") as fh:
        fh.write(rdir.name + "\n")
    link = root / "latest"
    try:
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(Path("releases") / rdir.name, link, target_is_directory=True)
    except OSError:  # e.g. Windows without symlink rights; LATEST file suffices
        pass


def latest_release(root: Path) -> str | None:
    f = root / "LATEST"
    return f.read_text().strip() if f.exists() else None
