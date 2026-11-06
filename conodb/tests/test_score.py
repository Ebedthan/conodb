"""Negatives and out-of-fold scoring."""

import random
import stat

import pytest

from conodb.config import load_config
from conodb.io import read_tsv
from conodb.negatives import shuffle_precursor
from conodb.pipeline import run


def test_shuffle_keeps_cysteines_and_region_composition():
    seq = "MKLTCVLIVAVLLLTACQLITAEDSRGTQKHRALRSGCCSDPRCAWRC"
    regions = [(1, 21), (22, 36), (37, 48)]
    out = shuffle_precursor(seq, regions, random.Random(1))
    assert out != seq and len(out) == len(seq)
    assert [i for i, c in enumerate(out) if c == "C"] == [i for i, c in enumerate(seq) if c == "C"]
    for s, e in regions:
        assert sorted(out[s - 1:e]) == sorted(seq[s - 1:e])
    assert shuffle_precursor(seq, regions, random.Random(1)) == out          # deterministic


def test_negatives(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    neg = {r["id"]: r for r in read_tsv(rdir / "negatives.tsv")}
    part = {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}
    assert neg["NGTOX1"]["role"] == "removed" and neg["NGTOX1"]["reason"].startswith(">= 0.9")
    assert neg["NGSHORT"]["role"] == "removed"
    real = [r for r in neg.values() if r["source"] == "uniprot" and r["role"] in ("train", "test")]
    assert len(real) == 24 and any(r["role"] == "test" for r in real)
    assert {r["source"] for r in neg.values()} >= {"uniprot", "venom_other", "shuffled"}
    # shuffled copies sit on their source precursor's side
    for r in neg.values():
        if r["source"] == "shuffled":
            src = part[r["accession"]]["role"]
            assert r["role"] == ("train" if src == "train" else "test")


def test_scores_are_out_of_fold(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    feats = {r["id"]: r for r in read_tsv(rdir / "features.tsv")}
    scores = read_tsv(rdir / "scores.tsv")
    part = {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}
    # every scored positive/negative is present, labels right
    assert {f["kind"] for f in feats.values()} == {"positive", "negative"}
    assert all(f["label"] == "NONE" for f in feats.values() if f["kind"] == "negative")
    # a training A entry: its own-family score comes from its fold's holdout
    a_train = sorted(a for a, p in part.items() if p["superfamily"] == "A" and p["role"] == "train")
    acc = a_train[0]
    lib = {s["library"] for s in scores if s["id"] == acc and s["superfamily"] == "A" and s["region"] == "SIG"}
    assert f"holdout:fold{part[acc]['fold']}" in lib
    oof = float(next(s["score"] for s in scores if s["id"] == acc and s["profile"] == "CONOPEP_A_SIG"
                     and s["library"].startswith("holdout")))
    assert float(feats[acc]["A:SIG:hmm"]) == oof
    # out-of-fold scores are lower on average than in-sample ones (no leakage)
    ins = {s["id"]: float(s["score"]) for s in scores if s["profile"] == "CONOPEP_A_SIG" and s["library"] == "final"}
    diffs = [ins[a] - float(feats[a]["A:SIG:hmm"]) for a in a_train if feats[a]["A:SIG:hmm"] not in ("", "NA")]
    assert diffs and sum(diffs) / len(diffs) > 0
    # tier-2 entries whose holdout could not be built are NA, not 0
    na = [f for f in feats.values() if "NA" in f.values()]
    assert all(f["role"] == "train" and f["fold"] == "logo" for f in na)
    # positives score higher on their own family's signal profile than negatives do
    pos = [float(feats[a]["A:SIG:hmm"]) for a in feats if feats[a]["label"] == "A" and feats[a]["A:SIG:hmm"] not in ("", "NA")]
    neg = [float(f["A:SIG:hmm"]) if f["A:SIG:hmm"] else 0.0 for f in feats.values() if f["kind"] == "negative"
           and f["source"] == "uniprot"]
    assert min(pos) > sorted(neg)[len(neg) // 2]


# The format pfsearch 2.3 prints with -f -k -a -r (as seen on a real run),
# alternating with the "/start-end" variant to cover both.
FAKE_PFSEARCH = r"""#!/bin/sh
for a in "$@"; do prf=$fa; fa=$a; done
name=$(awk '/^ID/{sub(";","",$2); print $2; exit}' "$prf")
awk -v base=100 -v m="$name" '/^>/{id=substr($1,2); n++;
  if (n % 2) printf(">%s motif=%s|%s raw_score=%d \n", id, m, m, base+n);
  else printf(">%s/1-10 level=NA raw_score=%d seq_end=-1\nXXXXXXXXXX\n", id, base+n)}' "$fa"
"""


def test_pssm_scoring_with_fake_pftools(project, mock_uniprot, fake_pftools):
    p = fake_pftools / "pfsearch"
    p.write_text(FAKE_PFSEARCH)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    toml = (project / "conodb.toml").read_text() + '\n[build]\npssm = "true"\n\n[score]\npssm = "true"\n'
    (project / "conodb.toml").write_text(toml)
    rdir = run(load_config(project / "conodb.toml"))
    feats = read_tsv(rdir / "features.tsv")
    assert "A:SIG:pssm" in feats[0]
    assert all(f["A:SIG:pssm"] not in ("",) for f in feats)
    pssm = [s for s in read_tsv(rdir / "scores.tsv") if s["method"] == "pssm"]
    assert pssm and all(int(s["score"]) > 100 for s in pssm)


def test_pfsearch_output_not_understood(project, mock_uniprot, fake_pftools):
    p = fake_pftools / "pfsearch"
    p.write_text("#!/bin/sh\necho 'something else'\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    toml = (project / "conodb.toml").read_text() + '\n[build]\npssm = "true"\n\n[score]\npssm = "true"\n'
    (project / "conodb.toml").write_text(toml)
    from conodb.score import ScoreError
    with pytest.raises(ScoreError, match="output format not understood"):
        run(load_config(project / "conodb.toml"))


def test_negative_cap(project, mock_uniprot):
    toml = (project / "conodb.toml").read_text() + "\n[negatives]\nmax_entries = 10\n"
    (project / "conodb.toml").write_text(toml)
    rdir = run(load_config(project / "conodb.toml"))
    neg = read_tsv(rdir / "negatives.tsv")
    capped = [r for r in neg if r["reason"].startswith("over negatives.max_entries")]
    kept = [r for r in neg if r["source"] == "uniprot" and r["role"] in ("train", "test")]
    assert len(capped) == 25 - 10 and len(kept) <= 10    # 26 minus one too short
