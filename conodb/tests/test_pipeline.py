"""End-to-end: fetch from a mock UniProt, parse, cluster, audit, reruns, updates."""

import json

import pytest

from conodb.cli import CURATION, OK, UPDATE_AVAILABLE, main
from conodb.config import load_config
from conodb.io import read_json, read_tsv
from conodb.pipeline import run
from conodb.uniprot import FetchError, fetch_all

from .conftest import MockUniProt, build_dataset, make_entry


def _families(rdir):
    """(superfamily, family) -> row; the conopeptide row wins over e.g. 'reference'."""
    rows = read_tsv(rdir / "families.tsv")
    out = {(r["superfamily"], r["family"]): r for r in rows if r["scope"] != "conopeptide"}
    out.update({(r["superfamily"], r["family"]): r for r in rows if r["scope"] == "conopeptide"})
    return out


def _ran(rdir):
    return {p.stem: read_json(p)["finished"] for p in (rdir / ".stamps").glob("*.json")}


def test_full_build(project, mock_uniprot, dataset):
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg)
    assert rdir.name == "2026_05"
    assert (project / "db" / "LATEST").read_text().strip() == "2026_05"

    entries = read_tsv(rdir / "entries.tsv")
    assert len(entries) == len(dataset)
    fams = _families(rdir)
    # A: 12 distinct + 1 near-duplicate (+1 conflicting copy labelled O1)
    assert fams[("A", "*")]["entries"] == "13"
    assert fams[("A", "*")]["tier"] == "1"
    assert int(fams[("A", "*")]["entry_groups"]) == 12
    # M: fragment shares the group of its precursor -> 3 groups -> tier 2
    assert fams[("M", "*")]["entries"] == "4" and fams[("M", "*")]["entry_groups"] == "3"
    assert fams[("M", "*")]["tier"] == "2"
    # nested label: C superfamily with Consomatin family
    assert fams[("C", "Consomatin")]["entries"] == "2"
    assert fams[("C", "Consomatin")]["tier"] == "3"   # family level: 2 groups
    assert fams[("C", "*")]["tier"] == "2"            # superfamily level: 3 groups
    # override maps insulin family into the classifier scope
    assert fams[("Con-insulin", "*")]["scope"] == "conopeptide"
    # uncurated labels are reported, not silently used
    unmapped = {r["raw_superfamily"] for r in read_tsv(rdir / "unmapped_labels.tsv")}
    assert unmapped == {"conoGAY family"}          # actinoporin is curated as 'other' in the fixture
    assert fams[("conoGAY", "*")]["tier"] == ""
    # QC: identical sequence labelled A and O1
    conflicts = read_tsv(rdir / "label_conflicts.tsv")
    assert any("PZ001:O1" in c["accessions"] and "PA001:A" in c["accessions"] for c in conflicts)
    # regions
    regions = read_tsv(rdir / "regions.tsv")
    pj = [r["region"] for r in regions if r["accession"] == "PJ001"]
    assert pj == ["SIG", "PRO", "MAT", "POST"]
    pq = [(r["region"], r["exact"]) for r in regions if r["accession"] == "PQ001"]
    assert pq == [("SIG", "true"), ("PRO", "false"), ("MAT", "true")]
    # short identical C-terminal propeptides are not reported as conflicts
    assert all("POST" not in c["kind"] for c in conflicts)
    fams_sens = _families(rdir)[("A", "*")]
    assert int(fams_sens["groups_0.9"]) >= int(fams_sens["groups_0.7"]) >= int(fams_sens["groups_0.5"])
    report = (rdir / "report.md").read_text()
    assert "Labels needing curation" in report and "conoGAY family" in report
    manifest = read_json(rdir / "manifest.json")
    assert manifest["uniprot"]["entries"] == len(dataset)
    assert set(manifest["stages"]) == {"fetch", "parse", "cluster", "audit", "partition", "negatives", "build", "score", "bayes", "fit",
                                        "select", "anatomy", "evaluate", "package"}
    # taxon scope: conotoxin-labelled turrid entry is a reference, not a class member
    assert {e["accession"]: e["scope"] for e in entries}["PT001"] == "reference"
    assert "Turridae" in report


def test_rerun_is_incremental(project, mock_uniprot):
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg)
    first = _ran(rdir)
    requests = mock_uniprot.requests
    run(cfg)
    assert _ran(rdir) == first                      # nothing reran
    assert mock_uniprot.requests == requests + 1    # only the release probe

    # changing tier thresholds reruns the audit only
    toml = (project / "conodb.toml").read_text().replace("tier1_min_groups = 10", "tier1_min_groups = 20")
    (project / "conodb.toml").write_text(toml)
    run(load_config(project / "conodb.toml"))
    after = _ran(rdir)
    assert sorted(s for s in after if after[s] != first[s]) == ["anatomy", "audit", "bayes", "build", "evaluate", "fit", "negatives", "package", "partition", "score", "select"]
    assert _families(rdir)[("A", "*")]["tier"] == "2"


def test_labels_file_change_reruns_parse(project, mock_uniprot):
    cfg = load_config(project / "conodb.toml")
    rdir = run(cfg)
    with open(project / "labels.tsv", "a") as fh:
        fh.write("conoGAY family\t*\tConoGAY\t\tconopeptide\t\n")
    run(load_config(project / "conodb.toml"))
    assert "conoGAY family" not in {r["raw_superfamily"] for r in read_tsv(rdir / "unmapped_labels.tsv")}
    assert _families(rdir)[("ConoGAY", "*")]["tier"] == "3"


def test_new_release_diff(project, dataset):
    m1 = MockUniProt(dataset, release="2026_04")
    _point(project, m1.url)
    run(load_config(project / "conodb.toml"))
    m1.close()

    newer = build_dataset()
    newer = [e for e in newer if e["primaryAccession"] != "PO004"]
    newer.append(make_entry("PN001", "Belongs to the conotoxin O1 superfamily.",
                            [("Signal", "MKLTCVLIVAVLLLTACQLIT", {}), ("Peptide", "DCCGVKLEMCHPCLC", {})]))
    m2 = MockUniProt(newer, release="2026_05")
    _point(project, m2.url)
    assert main(["check", "-c", str(project / "conodb.toml")]) == UPDATE_AVAILABLE
    rdir = run(load_config(project / "conodb.toml"))
    m2.close()
    assert rdir.name == "2026_05"
    diff = {(d["accession"], d["change"]) for d in read_tsv(rdir / "diff_entries.tsv")}
    assert diff == {("PN001", "added"), ("PO004", "removed")}
    assert read_json(rdir / "manifest.json")["previous_release"] == "2026_04"
    assert "Changes since 2026_04" in (rdir / "report.md").read_text()
    assert (project / "db" / "LATEST").read_text().strip() == "2026_05"


def test_strict_labels_stop_the_run(project, mock_uniprot):
    toml = (project / "conodb.toml").read_text().replace('overrides = "labels.tsv"', 'overrides = "labels.tsv"\nstrict = true')
    (project / "conodb.toml").write_text(toml)
    assert main(["run", "-c", str(project / "conodb.toml")]) == CURATION
    rdir = project / "db" / "releases" / "2026_05"
    assert (rdir / "unmapped_labels.tsv").exists()          # curation list is written
    assert not (project / "db" / "LATEST").exists()          # but nothing is published


def test_retries_and_truncation(tmp_path, dataset):
    cfg = load_config(None).source
    m = MockUniProt(dataset, fail_first=2)
    cfg.base_url, cfg.page_size, cfg.retries = m.url, 10, 4
    info = fetch_all(cfg, tmp_path / "x.jsonl.gz")
    assert info.total == len(dataset)
    m.close()
    m = MockUniProt(dataset, lie_total=len(dataset) + 5)
    cfg.base_url = m.url
    with pytest.raises(FetchError, match="incomplete"):
        fetch_all(cfg, tmp_path / "y.jsonl.gz")
    assert not (tmp_path / "y.jsonl.gz").exists()           # no partial file left
    m.close()


def test_offline_input(tmp_path, dataset):
    import gzip
    src = tmp_path / "dump.jsonl.gz"
    with gzip.open(src, "wt") as fh:
        for e in dataset:
            fh.write(json.dumps(e) + "\n")
    (tmp_path / "conodb.toml").write_text('[output]\nroot = "db"\n')
    assert main(["run", "-c", str(tmp_path / "conodb.toml"), "--input", str(src), "--release", "test"]) == OK
    assert (tmp_path / "db" / "releases" / "test" / "report.md").exists()


def _point(project, url):
    import re
    toml = (project / "conodb.toml").read_text()
    (project / "conodb.toml").write_text(re.sub(r'base_url = "[^"]+"', f'base_url = "{url}"', toml))


# ----------------------------------------------------------------- partition

def _partition(rdir):
    return {r["accession"]: r for r in read_tsv(rdir / "partition.tsv")}


def test_partition_design(project, mock_uniprot):
    rdir = run(load_config(project / "conodb.toml"))
    part = _partition(rdir)
    # tie conflict (A x1 vs O1 x1): both entries left out
    assert part["PA001"]["role"] == part["PZ001"]["role"] == "excluded"
    # A: 12 groups minus the excluded one -> tier 1, 20% test, 5 folds
    a = [r for r in part.values() if r["superfamily"] == "A" and r["scope"] == "conopeptide" and r["role"] != "excluded"]
    groups = {r["group"]: r["role"] for r in a}
    assert {r["tier"] for r in a} == {"1"}
    assert sum(v == "test" for v in groups.values()) == round(0.2 * len(groups))
    assert {r["fold"] for r in a if r["role"] == "train"} <= {"1", "2", "3", "4", "5"}
    # tier 2 -> leave-one-group-out, tier 3 -> reference, other scopes keep their role
    assert {r["fold"] for r in part.values() if r["superfamily"] == "O1" and r["role"] == "train"} == {"logo"}
    assert part["PJ001"]["role"] == "reference"
    assert part["PT001"]["role"] == "reference"
    assert part["PX000"]["role"] in ("unmapped", "negative")
    # no group on both sides
    sides = {}
    for r in part.values():
        if r["role"] in ("train", "test", "temporal_test"):
            sides.setdefault(r["group"], set()).add(r["role"] == "train")
    assert all(len(v) == 1 for v in sides.values())
    # deterministic
    before = (rdir / "partition.tsv").read_bytes()
    run(load_config(project / "conodb.toml"), force=True)
    assert (rdir / "partition.tsv").read_bytes() == before
    assert "Excluded entries (2)" in (rdir / "partition.md").read_text()


def test_partition_is_stable_across_releases(project, dataset):
    import random as _r
    from .conftest import _mutate, _rand

    m1 = MockUniProt(dataset, release="2026_04")
    _point(project, m1.url)
    old = _partition(run(load_config(project / "conodb.toml")))
    m1.close()

    rng = _r.Random(99)
    newer = [e for e in build_dataset() if e["primaryAccession"] != "PA005"]
    a0 = newer[0]
    sig = a0["sequence"]["value"][:21]
    for i in range(6):   # six new, distinct A groups
        newer.append(make_entry(f"PA{500 + i}", "Belongs to the conotoxin A superfamily.",
                                [("Signal", _mutate(rng, sig, 0.1), {}), ("Propeptide", _rand(rng, 24), {}),
                                 ("Peptide", _rand(rng, 0, "CC-C-C"), {})], first_public="2026-10-01"))
    m2 = MockUniProt(newer, release="2026_05")
    _point(project, m2.url)
    rdir = run(load_config(project / "conodb.toml"))
    m2.close()
    new = _partition(rdir)
    for acc, r in old.items():
        if acc in new and r["role"] in ("train", "test"):
            assert (new[acc]["role"], new[acc]["fold"]) == (r["role"], r["fold"]), acc
            assert new[acc]["inherited"] == "true"
    added = [new[f"PA{500 + i}"] for i in range(6)]
    assert all(r["inherited"] == "false" and r["role"] in ("train", "test") for r in added)
    assert "Previous release: 2026_04" in (rdir / "partition.md").read_text()


def test_temporal_test(project, mock_uniprot):
    toml = (project / "conodb.toml").read_text() + '\n[partition]\ntemporal_cutoff = "2020-01-01"\n'
    (project / "conodb.toml").write_text(toml)
    rdir = run(load_config(project / "conodb.toml"))
    part = _partition(rdir)
    # all fixture entries are dated 2012 except none -> no temporal test
    assert not any(r["role"] == "temporal_test" for r in part.values())
    toml = toml.replace('temporal_cutoff = "2020-01-01"', 'temporal_cutoff = "2010-01-01"')
    (project / "conodb.toml").write_text(toml)
    part = _partition(run(load_config(project / "conodb.toml")))
    roles = {r["role"] for r in part.values() if r["tier"] in ("1", "2")}
    assert roles == {"temporal_test"}
