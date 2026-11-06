"""Fit stage: design matrix, gradient, baseline, penalised and Stan fits."""


import numpy as np
import pytest
from scipy.optimize import check_grad

from conodb.cli import ERROR, main
from conodb.config import load_config
from conodb.fit import NONE, cmdstan_available, make_design, map_objective
from conodb.io import read_json, read_tsv
from conodb.pipeline import run


def test_design_semantics():
    cfg = load_config(None)
    fields = ["id", "kind", "source", "role", "label", "tier", "fold", "group",
              "A:SIG:hmm", "A:MAT:hmm", "B:SIG:hmm"]
    rows = [
        {"id": "a1", "kind": "positive", "source": "uniprot", "role": "train", "label": "A", "tier": "1",
         "fold": "1", "group": "g1", "A:SIG:hmm": "30", "A:MAT:hmm": "NA", "B:SIG:hmm": ""},
        {"id": "a2", "kind": "positive", "source": "uniprot", "role": "train", "label": "A", "tier": "1",
         "fold": "2", "group": "g1", "A:SIG:hmm": "20", "A:MAT:hmm": "5", "B:SIG:hmm": "1"},
        {"id": "n1", "kind": "negative", "source": "shuffled", "role": "test", "label": NONE, "tier": "",
         "fold": "", "group": "s", "A:SIG:hmm": "", "A:MAT:hmm": "", "B:SIG:hmm": ""},
    ]
    d = make_design(cfg, rows, fields)
    assert d.families == ["A", "B"] and d.feats == [("SIG", "hmm"), ("MAT", "hmm")]
    assert d.M[0, 0, 1] == 1 and np.isnan(d.raw[0, 0, 1])          # NA -> missing-profile indicator
    assert d.raw[0, 1, 0] == cfg.fit.hmm_floor                     # no score -> floor (-20 bits)
    assert np.isnan(d.raw[0, 1, 1]) and d.X[0, 1, 1] == 0           # B has no MAT profile: structural
    assert d.available.tolist() == [[True, True], [True, False]]
    assert list(d.y) == [0, 0, 2]
    assert d.w[0] == d.w[1] == 0.5 and d.w[2] == 0                  # weights only on train rows, 1/group size


def test_map_gradient():
    rng = np.random.default_rng(0)
    N, K, P = 30, 3, 2
    X, M = rng.normal(size=(N, K, P)), (rng.random((N, K, P)) < 0.2).astype(float)
    y, w = rng.integers(0, K + 1, N), rng.random(N)
    th = rng.normal(size=1 + K + P + K * P + P) * 0.3
    args = (X, M, y, w, K, P, 2.5, 0.5, 2.5)
    err = check_grad(lambda t: map_objective(t, *args)[0], lambda t: map_objective(t, *args)[1], th)
    assert err < 1e-4


def test_fit_stage(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    preds = read_tsv(rdir / "predictions.tsv")
    assert preds and {p["role"] for p in preds} <= {"test", "temporal_test"}
    model = read_json(rdir / "fit" / "model_map.json")
    assert model["classes"][-1] == NONE and model["converged"]
    # the fixture is easy: the penalised model gets the test set right
    acc = np.mean([p["map_pred"] == (p["label"] if p["kind"] == "positive" else NONE) for p in preds])
    assert acc >= 0.9
    # probabilities are probabilities
    assert all(0 <= float(p["map_prob"]) <= 1 for p in preds)
    cut = read_tsv(rdir / "cutoffs.tsv")
    assert cut and all(float(c["quantile"]) > 0.99 for c in cut)   # Bonferroni over profiles
    md = (rdir / "fit.md").read_text()
    assert "bayes: no posterior" in md and "Test set" in md


def test_bayes_required_without_cmdstan(project, mock_uniprot):
    toml = (project / "conodb.toml").read_text().replace('bayes = "false"', 'bayes = "true"\ncmdstan = "/nonexistent"')
    (project / "conodb.toml").write_text(toml)
    assert main(["run", "-c", str(project / "conodb.toml")]) == ERROR


@pytest.mark.skipif(not cmdstan_available(load_config(None)), reason="CmdStan not installed")
def test_bayes_fit(project, mock_uniprot):
    toml = (project / "conodb.toml").read_text().replace(
        'bayes = "false"', 'bayes = "true"\nchains = 2\niter_warmup = 300\niter_sampling = 300')
    (project / "conodb.toml").write_text(toml)
    rdir = run(load_config(project / "conodb.toml"))
    m = read_json(rdir / "fit" / "model_bayes.json")
    assert m["diagnostics"]["draws"] == 600
    assert (rdir / "fit" / "bayes_draws.npz").exists() and (rdir / "fit" / "bayes_summary.csv").exists()
    preds = read_tsv(rdir / "predictions.tsv")
    agree = np.mean([p["map_pred"] == p["bayes_pred"] for p in preds])
    assert agree >= 0.95, agree
    assert "Shared weights" in (rdir / "fit.md").read_text()


def test_negative_sources_balanced():
    cfg = load_config(None)
    fields = ["id", "kind", "source", "role", "label", "tier", "fold", "group", "A:SIG:hmm"]
    base = {"tier": "", "fold": "", "A:SIG:hmm": "1"}
    rows = [{"id": "p", "kind": "positive", "source": "uniprot", "role": "train", "label": "A", "group": "g"} | base]
    rows += [{"id": f"s{i}", "kind": "negative", "source": "shuffled", "role": "train", "label": NONE,
              "group": f"S{i}"} | base for i in range(9)]
    rows += [{"id": "r", "kind": "negative", "source": "uniprot", "role": "train", "label": NONE, "group": "R"} | base]
    d = make_design(cfg, rows, fields)
    w = dict(zip(d.ids, d.w, strict=True))
    assert abs(sum(w[f"s{i}"] for i in range(9)) - w["r"]) < 1e-9
    assert w["p"] == 1.0


def test_plan_cores():
    from conodb.fit import plan_cores
    assert plan_cores(16, 4) == (4, 4)
    assert plan_cores(8, 4) == (4, 2)
    assert plan_cores(6, 4) == (4, 1)
    assert plan_cores(2, 4) == (2, 1)
    assert plan_cores(1, 4) == (1, 1)
    assert plan_cores(12, 1) == (1, 12)


def test_evaluate_stage(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    md = (rdir / "evaluation.md").read_text()
    assert "By sequence type" in md and "Variant grid on the test set" in md and "Calibration" in md
    m = read_tsv(rdir / "evaluation_metrics.tsv")
    assert {r["group"] for r in m} >= {"all", "sequence type", "identity to training"}
    ab = read_tsv(rdir / "evaluation_ablations.tsv")
    names = {r["variant"] for r in ab}
    notes = {n for r in ab for n in r["note"].split("; ") if n}
    # the fixture has no pftools: PSSM variants are not run; the others are
    assert {"censored · HMM", "censored · HMM+PWM+FW", "floor -20 · HMM", "floor -20 · HMM+PWM"} <= names
    assert {"current settings", "signal only", "without signal"} <= notes
    assert not any("PSSM" in v for v in names)
    assert len(names) == len({(r["features"], r["hmm_floor"]) for r in ab})      # no duplicates


def test_design_floor_and_fill():
    from dataclasses import replace
    cfg = load_config(None)
    fields = ["id", "kind", "source", "role", "label", "tier", "fold", "group", "A:SIG:hmm", "A:SIG:pwm"]
    base = {"kind": "positive", "source": "uniprot", "role": "train", "label": "A", "tier": "1", "fold": "1"}
    rows = [base | {"id": "a", "group": "g1", "A:SIG:hmm": "-35", "A:SIG:pwm": "4.5"},
            base | {"id": "b", "group": "g2", "A:SIG:hmm": "-5", "A:SIG:pwm": ""},
            base | {"id": "c", "group": "g3", "A:SIG:hmm": "12", "A:SIG:pwm": "-2.0"}]
    cfg = replace(cfg, fit=replace(cfg.fit, methods=["hmm", "pwm"], hmm_floor=-20.0))
    d = make_design(cfg, rows, fields)
    assert d.raw[:, 0, 0].tolist() == [-20.0, -5.0, 12.0]           # clipped at the floor
    assert d.raw[1, 0, 1] == -2.0 and d.fill[1] == -2.0             # missing PWM -> lowest seen
    censored = make_design(replace(cfg, fit=replace(cfg.fit, hmm_floor=0.0)), rows, fields)
    assert censored.raw[:, 0, 0].tolist() == [0.0, 0.0, 12.0]


def test_apply_model_matches_design():
    """conodictor's core.apply_model gives the fit's probabilities for the same raw scores."""
    from conodb.core import apply_model
    from conodb.fit import export_model, predict_proba
    cfg = load_config(None)
    fields = ["id", "kind", "source", "role", "label", "tier", "fold", "group", "A:SIG:hmm", "B:SIG:hmm"]
    base = {"kind": "positive", "source": "uniprot", "role": "train", "tier": "1", "fold": "1"}
    rows = [base | {"id": f"a{i}", "group": f"g{i}", "label": "A", "A:SIG:hmm": str(10 + i), "B:SIG:hmm": ""} for i in range(3)]
    rows += [base | {"id": f"b{i}", "group": f"h{i}", "label": "B", "A:SIG:hmm": "-3", "B:SIG:hmm": str(15 + i)} for i in range(3)]
    d = make_design(cfg, rows, fields)
    rng = np.random.default_rng(0)
    alpha, B = rng.normal(size=2), rng.normal(size=(2, 1))
    model = export_model(d, alpha, B)
    p_fit = predict_proba(d, alpha, B)
    for i in range(len(d.ids)):
        raw = np.array([[float(rows[i]["A:SIG:hmm"]) if rows[i]["A:SIG:hmm"] else np.nan],
                        [float(rows[i]["B:SIG:hmm"]) if rows[i]["B:SIG:hmm"] else np.nan]])
        out = apply_model(model, raw)
        assert np.allclose(list(out["probabilities"].values()), p_fit[i], atol=1e-6)
