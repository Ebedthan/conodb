"""Build stage: leakage guards, usable profiles, reproducibility, PSSM plumbing."""


import pyhmmer
from pyhmmer import easel, plan7

from conodb.cli import ERROR, main
from conodb.config import load_config
from conodb.io import read_json, read_tsv
from conodb.pipeline import run


def _msa_accessions(path):
    return {line[1:].split("|")[0] for line in path.read_text().splitlines() if line.startswith(">")}


def test_profiles_never_contain_held_out_data(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    part = {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}
    rows = read_tsv(rdir / "profiles.tsv")
    built = [r for r in rows if r["status"] == "built"]
    assert built and (rdir / "profiles" / "final.hmm").exists()
    # A is tier 1 -> final + 5 fold holdouts; O1 is tier 2 -> one holdout per group
    a_units = {r["unit"] for r in built if r["superfamily"] == "A" and r["region"] == "SIG" and r["set"] != "ship"}
    assert a_units == {"", "fold1", "fold2", "fold3", "fold4", "fold5"}
    o1_units = {r["unit"] for r in built if r["superfamily"] == "O1" and r["region"] == "SIG" and r["set"] != "ship"}
    o1_groups = {p["group"] for p in part.values() if p["superfamily"] == "O1" and p["role"] == "train"}
    assert len(o1_units) == len(o1_groups) + 1
    for r in built:
        accs = _msa_accessions(rdir / r["msa"])
        if r["set"] == "ship":     # shipped profiles: train + test, never excluded/reference
            assert all(part[a]["role"] in ("train", "test", "temporal_test") for a in accs), r["name"]
            continue
        assert all(part[a]["role"] == "train" for a in accs), (r["name"], r["set"])
        if r["unit"].startswith("fold"):
            assert all(part[a]["fold"] != r["unit"][4:] for a in accs)
    # short C-terminal propeptides (< min_length) are not profiled
    assert not any(r["region"] == "POST" for r in built)
    # the shipped set includes held-out entries, the final set does not
    ship_a = _msa_accessions(rdir / next(r["msa"] for r in built if r["set"] == "ship" and r["name"] == "CONOPEP_A_SIG"))
    assert any(part[a]["role"] == "test" for a in ship_a)
    # full-precursor models: labels for every match state
    import json
    full = [r for r in built if r["region"] == "FULL"]
    assert full and all(len(json.load(open(rdir / r["labels"]))["labels"]) == int(r["match_length"]) for r in full)


def test_held_out_sequence_hits_its_own_family(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    part = {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}
    entries = {e["accession"]: e for e in read_tsv(rdir / "entries.tsv")}
    test_acc = sorted(a for a, p in part.items() if p["role"] == "test")
    assert test_acc
    alpha = easel.Alphabet.amino()
    with plan7.HMMFile(rdir / "profiles" / "final.hmm") as fh:
        hmms = list(fh)
    seqs = easel.TextSequenceBlock(
        [easel.TextSequence(name=a.encode(), sequence=entries[a]["sequence"]) for a in test_acc]).digitize(alpha)
    best = {}
    for hits in pyhmmer.hmmsearch(hmms, seqs, Z=1000, domZ=1000, T=0, domT=0):
        q = hits.query.name if isinstance(hits.query.name, str) else hits.query.name.decode()
        for h in hits:
            n = h.name if isinstance(h.name, str) else h.name.decode()
            if h.score > best.get(n, ("", -1e9))[1]:
                best[n] = (q, h.score)
    for a in test_acc:
        assert best[a][0].startswith(f"CONOPEP_{part[a]['superfamily']}_"), (a, best.get(a))


def test_build_is_reproducible_and_cached(project, mock_uniprot):
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg)
    first = (rdir / "profiles.tsv").read_bytes()
    final = (rdir / "profiles" / "final.hmm").read_bytes()
    run(cfg, from_stage="build")
    assert (rdir / "profiles.tsv").read_bytes() == first
    assert (rdir / "profiles" / "final.hmm").read_bytes() == final
    stats = read_json(rdir / ".stamps" / "build.json")["stats"]
    assert stats["reused"] == stats["profiles"]          # nothing rebuilt
    # a deleted profile file is detected and rebuilt
    victim = next(r for r in read_tsv(rdir / "profiles.tsv") if r["status"] == "built")
    (rdir / victim["hmm"]).unlink()
    run(cfg)
    assert (rdir / victim["hmm"]).exists()


def test_pssm_plumbing(project, mock_uniprot, fake_pftools):
    toml = (project / "conodb.toml").read_text() + '\n[build]\npssm = "true"\n'
    (project / "conodb.toml").write_text(toml)
    rdir = run(load_config(project / "conodb.toml"))
    built = [r for r in read_tsv(rdir / "profiles.tsv") if r["status"] == "built"]
    assert all(r["pssm"] for r in built if r["region"] != "FULL")   # full-precursor models: HMM only
    prf = (rdir / "profiles" / "final.prf").read_text()
    assert "ID   CONOPEP_A_SIG; MATRIX." in prf and "AC   CONOPEP_A_SIG;" in prf
    assert "DT   built by conodb." in prf and "SEQUENCE_PROFILE" not in prf
    assert "/tmp" not in prf and 'CC   /GENERATED_BY="conodb: pfw -m | pfmake -m";' in prf
    # byte-identical when rebuilt from scratch
    before = prf
    import shutil
    shutil.rmtree(rdir / "profiles")
    run(load_config(project / "conodb.toml"))
    assert (rdir / "profiles" / "final.prf").read_text() == before


def test_pssm_required_but_missing(project, mock_uniprot, monkeypatch):
    monkeypatch.setenv("PATH", "/nonexistent")
    toml = (project / "conodb.toml").read_text() + '\n[build]\npssm = "true"\n'
    (project / "conodb.toml").write_text(toml)
    assert main(["run", "-c", str(project / "conodb.toml")]) == ERROR
