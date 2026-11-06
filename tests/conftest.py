"""Synthetic UniProt-shaped entries and a local mock of the UniProt REST API."""

from __future__ import annotations

import json
import random
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import os
import stat

import pytest

AA = "ACDEFGHIKLMNPQRSTVWY"
NONC = AA.replace("C", "")
LINEAGE = ["Eukaryota", "Metazoa", "Mollusca", "Gastropoda", "Neogastropoda", "Conoidea", "Conidae", "Conus"]


def _feat(ftype, start, end, code="ECO:0000255", smod="EXACT", emod="EXACT", desc=""):
    return {
        "type": ftype,
        "location": {"start": {"value": start, "modifier": smod}, "end": {"value": end, "modifier": emod}},
        "description": desc,
        "evidences": [{"evidenceCode": code}],
    }


def make_entry(acc, similarity, parts, *, species="Conus testus", flag="Precursor", extra_features=(), first_public="2012-01-01",
               lineage=None):
    """parts: list of (feature_type or None, sequence, kwargs). Features get
    coordinates from their position in the concatenated sequence."""
    seq, feats, pos = "", [], 1
    for ftype, s, kw in parts:
        if ftype:
            feats.append(_feat(ftype, pos, pos + len(s) - 1, **kw))
        seq += s
        pos += len(s)
    entry = {
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "primaryAccession": acc,
        "uniProtkbId": f"{acc}_CONTE",
        "entryAudit": {"firstPublicDate": first_public, "lastAnnotationUpdateDate": "2025-06-18",
                       "entryVersion": 20, "sequenceVersion": 1},
        "organism": {"scientificName": species, "taxonId": 1000 + sum(map(ord, species)) % 1000, "lineage": lineage or LINEAGE},
        "proteinDescription": {"recommendedName": {"fullName": {"value": f"Conotoxin {acc}"}}, "flag": flag},
        "comments": [{"commentType": "SIMILARITY", "texts": [{"value": similarity}]}] if similarity else [],
        "features": feats + list(extra_features),
        "sequence": {"value": seq, "length": len(seq)},
    }
    return entry


def _mutate(rng, s, rate):
    return "".join(rng.choice(NONC) if (c != "C" and rng.random() < rate) else c for c in s)


def _rand(rng, n, cys_pattern=None):
    if cys_pattern:  # e.g. "CC-C-C" -> random spacers
        out = rng.choice(NONC)
        for ch in cys_pattern:
            out += "C" if ch == "C" else "".join(rng.choice(NONC) for _ in range(rng.randint(3, 6)))
        return out + rng.choice(NONC)
    return "".join(rng.choice(NONC) for _ in range(n or 12))


def build_dataset(seed=7):
    rng = random.Random(seed)
    entries = []

    def family(prefix, n, similarity, framework, species_pool=3):
        sig = "M" + _rand(rng, 20)
        pro = _rand(rng, 24)
        out = []
        for i in range(n):
            out.append(make_entry(
                f"{prefix}{i:03d}", similarity,
                [("Signal", _mutate(rng, sig, 0.1), {}),
                 ("Propeptide", _mutate(rng, pro, 0.6), {}),
                 ("Peptide", _rand(rng, 0, framework), {"code": "ECO:0000269"})],
                species=f"Conus sp{i % species_pool}",
            ))
        return out

    a = family("PA", 12, "Belongs to the conotoxin A superfamily.", "CC-C-C", 5)
    # near-duplicate of PA000 (same entry group)
    dup = json.loads(json.dumps(a[0]))
    dup["primaryAccession"] = "PA900"
    s = dup["sequence"]["value"]
    dup["sequence"]["value"] = s[:-2] + ("G" if s[-2] != "G" else "A") + s[-1]
    a.append(dup)
    entries += a
    entries += family("PO", 5, "Belongs to the conotoxin O1 superfamily.", "C-C-CC-C-C")
    m = family("PM", 3, "Belongs to the conotoxin M superfamily.", "CC-C-C-CC")
    entries += m
    # mature-only fragment identical to PM000's mature peptide
    mat = [f for f in m[0]["features"] if f["type"] == "Peptide"][0]
    mseq = m[0]["sequence"]["value"][mat["location"]["start"]["value"] - 1: mat["location"]["end"]["value"]]
    entries.append(make_entry("PM900", "Belongs to the conotoxin M superfamily.",
                              [("Peptide", mseq, {"code": "ECO:0000269"})], flag="Fragment"))
    entries += family("PS", 2, "Belongs to the conotoxin C superfamily. Consomatin family.", "C-C")
    entries += family("PC", 1, "Belongs to the conotoxin C superfamily.", "")
    entries += family("PG", 1, "Belongs to the conoGAY family.", "C-C")
    entries += family("PI", 2, "Belongs to the insulin family.", "C-C-C-C")
    entries += family("PX", 1, "Belongs to the actinoporin family. Conoidea subfamily.", "")
    entries += family("PU", 1, "", "C-C")  # no SIMILARITY comment
    # uncertain propeptide boundary
    entries.append(make_entry("PQ001", "Belongs to the conotoxin O1 superfamily.",
                              [("Signal", "M" + _rand(rng, 20), {}),
                               ("Propeptide", _rand(rng, 20), {"emod": "UNSURE"}),
                               ("Peptide", _rand(rng, 0, "C-C-CC-C-C"), {})]))
    # chain-only mature + C-terminal propeptide (POST)
    entries.append(make_entry("PJ001", "Belongs to the conotoxin J superfamily.",
                              [("Signal", "M" + _rand(rng, 20), {}),
                               ("Propeptide", _rand(rng, 18), {}),
                               ("Chain", _rand(rng, 0, "C-C-C-C"), {}),
                               ("Propeptide", "GKR", {})]))
    # a conotoxin-labelled entry from outside Conidae (taxon scope -> reference)
    entries.append(make_entry("PT001", "Belongs to the conotoxin A superfamily.",
                              [("Signal", "M" + _rand(rng, 20), {}), ("Peptide", _rand(rng, 0, "CC-C-C"), {})],
                              species="Turris testus", lineage=LINEAGE[:-2] + ["Turridae", "Turris"]))
    # label conflict: same sequence as PA001 but labelled O1
    conflict = json.loads(json.dumps(a[1]))
    conflict["primaryAccession"] = "PZ001"
    conflict["comments"][0]["texts"][0]["value"] = "Belongs to the conotoxin O1 superfamily."
    entries.append(conflict)
    return entries


class MockUniProt:
    """Serves /uniprotkb/search with cursor pagination and UniProt headers."""

    def __init__(self, entries, release="2026_05", fail_first=0, lie_total=None, routes=None):
        self.entries = entries
        self.routes = routes or {}   # query substring -> entries (e.g. the negatives query)
        self.release = release
        self.fail_first = fail_first
        self.lie_total = lie_total
        self.requests = 0
        mock = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def do_GET(self):
                mock.requests += 1
                if mock.fail_first > 0:
                    mock.fail_first -= 1
                    self.send_response(503)
                    self.send_header("Retry-After", "0")
                    self.end_headers()
                    return
                u = urllib.parse.urlparse(self.path)
                q = urllib.parse.parse_qs(u.query)
                size = int(q.get("size", ["25"])[0])
                cursor = int(q.get("cursor", ["0"])[0])
                query = q.get("query", [""])[0]
                pool = next((v for k, v in mock.routes.items() if k in query), mock.entries)
                page = pool[cursor: cursor + size]
                body = json.dumps({"results": page}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("X-Total-Results", str(mock.lie_total or len(pool)))
                self.send_header("X-UniProt-Release", mock.release)
                self.send_header("X-UniProt-Release-Date", "05-November-2026")
                if cursor + size < len(pool):
                    params = dict((k, v[0]) for k, v in q.items())
                    params["cursor"] = str(cursor + size)
                    nxt = f"http://127.0.0.1:{self.server.server_port}{u.path}?{urllib.parse.urlencode(params)}"
                    self.send_header("Link", f'<{nxt}>; rel="next"')
                self.end_headers()
                self.wfile.write(body)

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.server.server_port}/uniprotkb/search"
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def close(self):
        self.server.shutdown()
        self.server.server_close()


@pytest.fixture
def dataset():
    return build_dataset()


def build_negatives_dataset(dataset, seed=11):
    """Secreted non-toxin proteins (signal peptide + body), plus two QC cases."""
    rng = random.Random(seed)
    out = []
    for i in range(24):
        body = "".join(rng.choice(AA) for _ in range(rng.randint(60, 180)))
        species, lineage = (("Conus sp", LINEAGE) if i % 4 == 0 else
                            ("Lottia gigantea", LINEAGE[:4] + ["Patellogastropoda", "Lottiidae", "Lottia"]))
        out.append(make_entry(f"NG{i:03d}", "", [("Signal", "M" + _rand(rng, 20), {}), (None, body, {})],
                              species=species, lineage=lineage))
    # identical to a conopeptide precursor -> removed by QC
    toxin_like = json.loads(json.dumps(dataset[3]))
    toxin_like["primaryAccession"] = "NGTOX1"
    toxin_like["comments"] = []
    out.append(toxin_like)
    out.append(make_entry("NGSHORT", "", [(None, "MKLLV", {})]))           # too short
    return out


@pytest.fixture
def negatives_dataset(dataset):
    return build_negatives_dataset(dataset)


@pytest.fixture
def mock_uniprot(dataset, negatives_dataset):
    m = MockUniProt(dataset, routes={"33208": negatives_dataset})
    yield m
    m.close()


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
    import conodb.uniprot as up
    monkeypatch.setattr(up.time, "sleep", lambda s: None)


LABELS_TSV = (
    "raw_superfamily\traw_family\tsuperfamily\tfamily\tscope\tnote\n"
    "# comment line\n"
    "insulin family\t*\tCon-insulin\t\tconopeptide\t\n"
    "actinoporin family\t*\tActinoporin\t\tother\t\n"
)


@pytest.fixture
def project(tmp_path, mock_uniprot):
    """A config dir pointing at the mock server."""
    (tmp_path / "labels.tsv").write_text(LABELS_TSV)
    (tmp_path / "conodb.toml").write_text(
        f'''[source]
query = "taxonomy_id:37797 AND reviewed:true"
base_url = "{mock_uniprot.url}"
page_size = 7
retries = 3

[labels]
overrides = "labels.tsv"
taxon_scope = ["Conidae"]

[cluster]
identities = [0.9, 0.7, 0.5]
group_identity = 0.7

[tiers]
tier1_min_groups = 10
tier2_min_groups = 3

[fit]
bayes = "false"

[output]
root = "db"
'''
    )
    return tmp_path


FAKE_PFW = """#!/bin/sh
[ "$1" = "-m" ] || { echo "expected -m" >&2; exit 2; }
sed 's/^>\\(.*\\)$/>\\1 weight= 1.000/' "$2"
"""
FAKE_PFMAKE = """#!/bin/sh
[ "$1" = "-m" ] || { echo "expected -m" >&2; exit 2; }
[ -f "$3" ] || { echo "no matrix $3" >&2; exit 2; }
grep -q 'weight=' "$2" || { echo "unweighted input" >&2; exit 2; }
cat <<PRF
ID   SEQUENCE_PROFILE; MATRIX.
AC   ?;
DT   05-NOV-2026
DE   Generated from MSA file.
MA   /GENERAL_SPEC: ALPHABET='ABCDEFGHIKLMNPQRSTVWYZ'; LENGTH=4;
CC   /GENERATED_BY="/usr/local/bin/pfmake -m $2 $3";
//
PRF
"""


@pytest.fixture
def fake_pftools(tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name, body in (("pfw", FAKE_PFW), ("pfmake", FAKE_PFMAKE)):
        p = bindir / name
        p.write_text(body)
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    share = tmp_path / "share" / "pftools"
    share.mkdir(parents=True)
    (share / "blosum45.cmp").write_text("matrix\n")
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("PFTOOLS_DATA", str(tmp_path / "share"))
    return bindir


