"""The Bayesian fit is expensive: it must only run when its posterior would change."""

import numpy as np
import pandas as pd
import pytest

import conodb.fit as fitmod
from conodb.config import load_config
from conodb.io import read_json, sha256_file, sha256_obj, write_json
from conodb.pipeline import run


def _enable_bayes(project, extra=""):
    toml = (project / "conodb.toml").read_text().replace('bayes = "false"', 'bayes = "true"' + extra)
    (project / "conodb.toml").write_text(toml)


@pytest.fixture
def fake_stan(monkeypatch):
    """Replace sampling by a quick fake that records every call."""
    calls = []

    def fake_fit_bayes(cfg, d, workdir, cache):
        calls.append(cfg.fit.chains)
        m = fitmod.fit_map(cfg, d)
        draws = 40
        alpha = np.repeat(m["alpha"][None], draws, 0)
        B = np.repeat(m["B"][None], draws, 0)
        return {"alpha_draws": alpha, "B_draws": B, "summary": pd.DataFrame({"R_hat": [1.0]}, index=["beta[1]"]),
                "max_rhat": 1.0, "min_ess_bulk": 500.0, "divergences": 0, "draws": draws,
                "beta_mean": m["beta"], "tau_mean": np.full(len(m["beta"]), 0.4),
                "parallel_chains": 1, "threads_per_chain": 1, "seconds": 0.1}

    monkeypatch.setattr(fitmod, "fit_bayes", fake_fit_bayes)
    monkeypatch.setattr(fitmod, "cmdstan_available", lambda cfg: True)
    import conodb.pipeline as pl
    monkeypatch.setattr(pl, "_cmdstan_id", lambda cfg: "/opt/cmdstan")
    return calls


def test_bayes_only_reruns_when_the_posterior_would_change(project, mock_uniprot, fake_stan):
    _enable_bayes(project)
    rdir = run(load_config(project / "conodb.toml"))
    assert fake_stan == [4]
    assert read_json(rdir / "fit" / "bayes_status.json")["status"] == "sampled"
    assert "bayes_pred" in (rdir / "predictions.tsv").read_text().splitlines()[0]

    # cores, grainsize, the fast fit's settings, a forced fit/evaluate rerun: no sampling
    toml = (project / "conodb.toml").read_text().replace("[fit]", "[fit]\ncores = 3\ngrainsize = 7\ntau_fixed = 0.3\nbaseline_fpr = 0.05")
    (project / "conodb.toml").write_text(toml)
    run(load_config(project / "conodb.toml"))
    run(load_config(project / "conodb.toml"), from_stage="fit")
    assert fake_stan == [4]
    # the stage reran (its key changed) but found the posterior and reused it
    run(load_config(project / "conodb.toml"), from_stage="bayes")
    assert fake_stan == [4]
    assert read_json(rdir / "fit" / "bayes_status.json")["status"] == "reused"

    # a sampler setting changes the posterior: sample again
    toml = (project / "conodb.toml").read_text().replace("[fit]", "[fit]\nchains = 2")
    (project / "conodb.toml").write_text(toml)
    run(load_config(project / "conodb.toml"))
    assert fake_stan == [4, 2]


def _old_fit_stage_key(cfg, rdir):
    """The key conodb 0.1.0's single 'fit' stage stamped (copied from that version)."""
    return sha256_obj({
        "stage": "fit",
        "conodb": "0.1.0",
        "config": {s: cfg.section(s) for s in ("fit", "score")},
        "inputs": {"features.tsv": sha256_file(rdir / "features.tsv")},
        "extra": [cfg.fit.bayes != "false" and "/opt/cmdstan", 4, sha256_file(fitmod.STAN_FILE)],
    })


def test_posterior_from_previous_version_is_adopted(project, mock_uniprot, fake_stan, monkeypatch):
    _enable_bayes(project)
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg, until_stage="score")
    # what the previous version left behind: draws, model_bayes.json without a key, and its fit stamp
    d = fitmod.load_design(cfg, rdir)
    K, P = len(d.families), len(d.feats)
    rng = np.random.default_rng(1)
    (rdir / "fit").mkdir(exist_ok=True)
    np.savez_compressed(rdir / "fit" / "bayes_draws.npz", alpha=rng.normal(size=(30, K)), B=rng.normal(size=(30, K, P)))
    write_json(rdir / "fit" / "model_bayes.json", {
        "families": d.families, "classes": d.families + ["NONE"], "features": [f"{r}:{m}" for r, m in d.feats],
        "beta_mean": [0.0] * P, "tau_mean": [0.5] * P,
        "diagnostics": {"max_rhat": 1.006, "min_ess_bulk": 796, "divergences": 0, "draws": 4000}})
    write_json(rdir / ".stamps" / "fit.json", {"key": _old_fit_stage_key(cfg, rdir), "outputs": {}})

    def no_sampling(*a, **k):
        raise AssertionError("must not sample")
    monkeypatch.setattr(fitmod, "fit_bayes", no_sampling)
    run(cfg)
    status = read_json(rdir / "fit" / "bayes_status.json")
    assert status["status"] == "adopted"
    assert "bayes" in (rdir / "fit.md").read_text() and "max R-hat 1.006" in (rdir / "fit.md").read_text()
    # and from then on it is simply reused
    run(cfg, from_stage="bayes")
    assert read_json(rdir / "fit" / "bayes_status.json")["status"] == "reused"


def test_stale_posterior_is_not_adopted(project, mock_uniprot, fake_stan):
    _enable_bayes(project)
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg, until_stage="score")
    (rdir / "fit").mkdir(exist_ok=True)
    np.savez_compressed(rdir / "fit" / "bayes_draws.npz", alpha=np.zeros((2, 1)), B=np.zeros((2, 1, 1)))
    write_json(rdir / "fit" / "model_bayes.json", {"families": ["X"], "features": ["SIG:hmm"]})
    write_json(rdir / ".stamps" / "fit.json", {"key": "something else", "outputs": {}})
    run(cfg)
    assert read_json(rdir / "fit" / "bayes_status.json")["status"] == "sampled"
    assert fake_stan == [4]
