"""Region prediction benchmark, framework feature, package, train/serve consistency."""

import json
import tarfile

import numpy as np

from conodb.config import load_config
from conodb.core import apply_model
from conodb.fit import bayes_core_key, make_design
from conodb.io import read_json, read_tsv
from conodb.pipeline import run


def test_anatomy_benchmark(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    seq = read_tsv(rdir / "anatomy_sequences.tsv")
    allrows = {(r["family_source"], r["variant"]): r for r in seq if r["subset"] == "all"}
    from conodb.variants import region_variants
    variants = {v["name"] for v in region_variants(load_config(project / "conodb.toml"))}
    assert len(variants) == 44
    assert {v for s, v in allrows if s == "true"} == variants
    assert {v for s, v in allrows if s == "predicted"} == variants
    # fixture precursors are clean: label transfer recovers most mature peptides exactly
    assert float(allrows[("true", "transfer/within/oh5")]["mature_exact"]) >= 0.8
    assert float(allrows[("true", "transfer/cross/oh5")]["mature_exact"]) >= 0.8
    assert (rdir / "anatomy_errors.tsv").exists()
    m = read_tsv(rdir / "anatomy_metrics.tsv")
    assert {r["region"] for r in m} == {"SIG", "PRO", "MAT", "POST"}       # POST evaluated too
    assert "Region prediction benchmark" in (rdir / "anatomy.md").read_text()


def test_framework_feature_is_out_of_fold(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    part = {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}
    fw = [s for s in read_tsv(rdir / "scores.tsv") if s["method"] == "framework"]
    assert fw
    for s in fw:
        p = part.get(s["id"])
        own = p and p["role"] == "train" and s["superfamily"] == p["superfamily"]
        assert s["library"].startswith("holdout:") == bool(own), s
    header = (rdir / "features.tsv").read_text().splitlines()[0].split("\t")
    assert "A:FW:framework" in header and "A:SIG:pwm" in header


def test_bayes_key_ignores_unused_columns():
    from dataclasses import replace
    cfg = load_config(None)
    cfg = replace(cfg, fit=replace(cfg.fit, methods=["hmm"]))
    fields = ["id", "kind", "source", "role", "label", "tier", "fold", "group", "A:SIG:hmm", "A:SIG:pwm"]
    base = {"kind": "positive", "source": "uniprot", "role": "train", "label": "A", "tier": "1", "fold": "1"}
    rows = [base | {"id": f"a{i}", "group": f"g{i}", "A:SIG:hmm": str(10 + i), "A:SIG:pwm": "1"} for i in range(3)]
    k1 = bayes_core_key(cfg, make_design(cfg, rows, fields))
    rows2 = [r | {"A:SIG:pwm": "7"} for r in rows]
    assert bayes_core_key(cfg, make_design(cfg, rows2, fields)) == k1          # unused column changed
    rows3 = [rows[0] | {"A:SIG:hmm": "30"}] + rows[1:]
    assert bayes_core_key(cfg, make_design(cfg, rows3, fields)) != k1          # used column changed


def test_package_and_consistency(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    info = read_json(rdir / "package.json")
    pkg = rdir / info["package"]
    for f in ("profiles.hmm", "precursors.hmm", "precursors.labels.json", "sites.json", "framework.json",
              "model.json", "settings.json", "manifest.json", "README.md", "cutoffs.tsv", "supplementary/README.md"):
        assert (pkg / f).exists(), f
    hmm_text = (pkg / "profiles.hmm").read_text()
    assert hmm_text.count("\nGA ") == hmm_text.count("\nNAME ")                  # every profile has GA
    with tarfile.open(rdir / info["archive"]) as tf:
        assert any(n.endswith("model.json") for n in tf.getnames())
    sup = (pkg / "supplementary" / "README.md").read_text()
    assert "evaluation_ablations.tsv" in sup and "anatomy_metrics.tsv" in sup and "conodb.toml" in sup
    manifest = read_json(pkg / "manifest.json")
    assert manifest["model"].startswith("conditional logit, penalised, refitted on all")

    # core.apply_model on raw features reproduces the evaluation model's predictions
    model = read_json(rdir / "fit" / "model_map.json")
    feats = model["features"]
    rows = {r["id"]: r for r in read_tsv(rdir / "features.tsv")}
    preds = read_tsv(rdir / "predictions.tsv")
    for p in preds[:25]:
        r = rows[p["id"]]
        raw = np.full((len(model["families"]), len(feats)), np.nan)
        for k, fam in enumerate(model["families"]):
            for j, f in enumerate(feats):
                v = r.get(f"{fam}:{f}")
                if v not in (None, "", "NA"):
                    raw[k, j] = float(v)
        out = apply_model(model, raw)
        assert out["prediction"] == p["map_pred"]
        assert abs(out["p_none"] - float(p["map_p_none"])) < 1e-3
        assert out["best_family"] == p["map_best_family"]
    json.dumps(model)


def test_cterm_trim():
    from conodb.core import trim_cterm
    seq = "MKLCCKLGRR"
    assert trim_cterm(seq, {"MAT": (1, 10)}, 3) == ({"MAT": (1, 7)}, True)          # GRR removed
    assert trim_cterm("MKLCCKLG", {"MAT": (1, 8)}, 3) == ({"MAT": (1, 7)}, True)
    assert trim_cterm("MKLCCKLA", {"MAT": (1, 8)}, 3)[1] is False
    assert trim_cterm(seq, {"MAT": (1, 10), "POST": (11, 14)}, 3)[1] is False       # something follows
    assert trim_cterm("GRR", {"MAT": (1, 3)}, 3)[1] is False                         # would leave nothing


def test_transfer_terminal_overhang():
    from conodb.core import DomainHit, transfer_labels
    lab = list("SSSPPMMMM")
    h = DomainHit(score=1, start=3, end=11, hmm_from=1, hmm_to=9, hmm_aln="x" * 9, seq_aln="A" * 9)
    assert "".join(transfer_labels(14, h, lab, 2)) == "SSSSSPPMMMMMM-"
    assert "".join(transfer_labels(14, h, lab, 0)) == "--SSSPPMMMM---"


def test_supplementary_names_unique(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    sup = rdir / read_json(rdir / "package.json")["package"] / "supplementary"
    nums = [f.name.split("_", 1)[0] for f in sup.iterdir() if f.name != "README.md"]
    assert len(nums) == len(set(nums))
