"""Selection stage: folds, rule, outputs."""

import numpy as np

from conodb.config import load_config
from conodb.io import read_json, read_tsv
from conodb.pipeline import run
from conodb.select import choose, cv_folds


def test_cv_folds_keep_groups_and_shuffled_with_source():
    rows = [
        {"id": "p1", "kind": "positive", "fold": "2", "group": "G1"},
        {"id": "p2", "kind": "positive", "fold": "2", "group": "G1"},
        {"id": "p3", "kind": "positive", "fold": "logo", "group": "G7"},
        {"id": "p4", "kind": "positive", "fold": "logo", "group": "G7"},
        {"id": "s1", "kind": "negative", "fold": "", "group": "SHUF-G1"},
        {"id": "n1", "kind": "negative", "fold": "", "group": "N3"},
        {"id": "n2", "kind": "negative", "fold": "", "group": "N3"},
    ]
    f = cv_folds(rows, 5)
    assert f["p1"] == f["p2"] == 2 == f["s1"]
    assert f["p3"] == f["p4"] and f["n1"] == f["n2"]
    assert all(1 <= v <= 5 for v in f.values())
    assert cv_folds(rows, 5) == f                       # stable


def test_choose_simplest_not_worse():
    point = np.array([0.30, 0.31, 0.50])
    rng = np.random.default_rng(0)
    base = rng.normal(0, 0.02, 500)
    boot = np.column_stack([0.30 + base, 0.31 + base + rng.normal(0, 0.02, 500), 0.50 + base])
    best, sel, diff, nw = choose(point, boot, [5, 2, 1], higher_better=False)
    assert best == 0 and sel == 1                        # 1 is not clearly worse and simpler; 2 is worse
    assert nw.tolist() == [True, True, False]


def test_choose_guard_protects_subset():
    point = np.array([0.30, 0.30])
    rng = np.random.default_rng(1)
    e = rng.normal(0, 0.005, 400)
    boot = np.column_stack([0.30 + e, 0.30 + e])
    gp = np.array([0.4, 1.5])                          # variant 1 (simpler) fails mature peptides
    gb = np.column_stack([0.4 + e, 1.5 + e])
    _, sel, _, elig = choose(point, boot, [3, 1], higher_better=False, guards=((gp, gb),))
    assert sel == 0 and elig.tolist() == [True, False]


def test_select_stage_outputs(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    rec = read_json(rdir / "selection_recommended.json")
    assert rec["fit"]["methods"] and "hmm_floor" in rec["fit"]
    assert rec["anatomy"]["method"] in ("transfer", "envelope")
    cls = read_tsv(rdir / "selection_classification.tsv")
    taus = {float(r["tau"]) for r in cls}
    assert taus == {0.25, 0.5, 1.0}
    assert sum(r["selected"] == "true" for r in cls if float(r["tau"]) == 0.5) == 1
    if any(r["methods"] == "hmm+pssm" for r in cls):
        assert any(r["note"].startswith("previous ConoDictor") for r in cls)
    reg = read_tsv(rdir / "selection_regions.tsv")
    assert sum(r["selected"] == "true" for r in reg) == 1
    md = (rdir / "selection.md").read_text()
    assert "Shrinkage sensitivity" in md and "Region prediction" in md
    ab = read_tsv(rdir / "evaluation_ablations.tsv")
    assert sum(r["recommended"] == "true" for r in ab if r["subset"] == "all") == 1
    assert all(r["log_loss_lo"] != "" for r in ab if r["subset"] == "all")
