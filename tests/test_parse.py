from conodb.config import LabelsConfig, RegionsConfig
from conodb.parse import LabelResolver, cysteine_framework, evidence_class, extract_regions, split_similarity

from .conftest import make_entry


def test_split_similarity():
    assert split_similarity("Belongs to the conotoxin A superfamily.") == ("conotoxin A superfamily", "")
    assert split_similarity("Belongs to the conotoxin C superfamily. Consomatin family.") == (
        "conotoxin C superfamily", "Consomatin family")
    assert split_similarity("Belongs to the actinoporin family. Conoidea subfamily.") == (
        "actinoporin family", "Conoidea subfamily")


def test_cysteine_framework():
    assert cysteine_framework("GCCSDPRCAWRC") == "CC-C-C"
    assert cysteine_framework("ZKKGCCSNPACRTNHPEVCX") == "CC-C-C"
    assert cysteine_framework("GEEHSKYQ") == ""
    assert cysteine_framework("dcCsk") == "CC"


def test_evidence_class():
    assert evidence_class(["ECO:0000255", "ECO:0000269"]) == "experimental"
    assert evidence_class(["ECO:0000250"]) == "similarity"
    assert evidence_class([]) == "none"


def test_label_rules_and_overrides(tmp_path):
    f = tmp_path / "labels.tsv"
    f.write_text("raw_superfamily\traw_family\tsuperfamily\tfamily\tscope\n"
                 "insulin family\t*\tCon-insulin\t\tconopeptide\n"
                 "conotoxin O1 superfamily\t*\tO1\t\texclude\n")
    r = LabelResolver.from_config(LabelsConfig(overrides=f))
    lab = r.resolve("conotoxin C superfamily", "Consomatin family")
    assert (lab.superfamily, lab.family, lab.scope, lab.source) == ("C", "Consomatin", "conopeptide", "rule")
    assert r.resolve("conopeptide P-like superfamily", "").superfamily == "P-like"
    assert r.resolve("conotoxin B2 family", "").superfamily == "B2"
    # bare codes, as UniProt writes some entries
    lab = r.resolve("O2 superfamily", "Contryphan family")
    assert (lab.superfamily, lab.family, lab.scope) == ("O2", "Contryphan", "conopeptide")
    assert r.resolve("E superfamily", "").superfamily == "E"
    assert r.resolve("Pg turripeptide superfamily", "").scope == "unmapped"
    assert r.resolve("teretoxin M (TM) superfamily", "").scope == "unmapped"
    assert r.resolve("insulin family", "").superfamily == "Con-insulin"
    assert r.resolve("conotoxin O1 superfamily", "").scope == "exclude"   # override beats rule
    lab = r.resolve("conoGAY family", "")
    assert lab.scope == "unmapped" and lab.superfamily == "conoGAY"
    assert r.unmapped[("conoGAY family", "")] == 1
    assert r.resolve("", "").scope == "unlabelled"


def _regions(entry, **kw):
    rows, warnings = extract_regions(entry, RegionsConfig(**kw))
    return [(r["region"], r["index"], r["sequence"]) for r in rows], warnings


def test_regions_pro_mat_post_inter():
    e = make_entry("X1", "", [
        ("Signal", "MKLTCVLIVAVLLLTACQLIT", {}),
        ("Propeptide", "AEDSRGTQKHRALRS", {}),
        ("Peptide", "GCCSDPRCAWRC", {}),
        ("Propeptide", "GRR", {}),
        ("Peptide", "GCCSDPRCNMNNPDYC", {}),
        ("Propeptide", "GKR", {}),
    ])
    rows, _ = _regions(e, min_length=3)
    assert [(r, i) for r, i, _ in rows] == [("SIG", 1), ("PRO", 1), ("MAT", 1), ("INTER", 1), ("MAT", 2), ("POST", 1)]


def test_chain_fallback_and_uncertain():
    e = make_entry("X2", "", [
        ("Signal", "MKLTCVLIVAVLLLTACQLIT", {}),
        ("Propeptide", "AEDSRGTQKHRALRS", {"emod": "UNSURE"}),
        ("Chain", "GCCSDPRCAWRC", {}),
    ])
    rows, warnings = _regions(e, allow_uncertain=False)
    assert [r for r, _, _ in rows] == ["SIG", "MAT"]
    assert any("uncertain boundary" in w for w in warnings)
    rows, warnings = _regions(e)  # default: kept, flagged exact=false
    assert [r for r, _, _ in rows] == ["SIG", "PRO", "MAT"] and not warnings
    full, _ = extract_regions(e, RegionsConfig())
    assert [r["exact"] for r in full] == [True, False, True]


def test_chain_spanning_precursor_is_not_mature():
    e = make_entry("X3", "", [("Signal", "MKLTCVLIVAVLLLTACQLIT", {}), (None, "AEDSRGTQKHRALRSGCCSDPRCAWRC", {})])
    e["features"].append({"type": "Chain", "location": {"start": {"value": 1, "modifier": "EXACT"},
                          "end": {"value": 48, "modifier": "EXACT"}}, "evidences": []})
    rows, _ = _regions(e)
    assert [r for r, _, _ in rows] == ["SIG"]


def test_entry_overrides(tmp_path):
    f = tmp_path / "entries.tsv"
    f.write_text("accession\tsuperfamily\tfamily\tscope\tnote\nP0DKQ3\tO1\t\tconopeptide\tsignal is O1\n"
                 "Q00001\tConoGAY\t\tconopeptide\t\n")
    r = LabelResolver.from_config(LabelsConfig(entry_overrides=f))
    lab = r.apply_entry("P0DKQ3", r.resolve("conotoxin M superfamily", ""))
    assert (lab.superfamily, lab.source) == ("O1", "entry_override")
    lab = r.apply_entry("Q00001", r.resolve("conoGAY family", ""))
    assert lab.superfamily == "ConoGAY" and r.unmapped == {}      # curated -> no longer unmapped
    assert r.apply_entry("Q99999", r.resolve("conotoxin A superfamily", "")).source == "rule"
