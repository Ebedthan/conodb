"""Shared core: everything conodictor needs at prediction time.

conodb trains with these functions and conodictor predicts with the very
same ones, so features are computed identically in both (no train/serve
skew). Nothing here depends on the pipeline, files of a release, or config
objects: inputs are sequences, models and plain parameters.

Coordinates are 1-based and inclusive, as in UniProt. A "segments" dict maps
a region (SIG, PRO, MAT, POST) to (start, end).

Region prediction
  transfer   align the sequence to a family's full-precursor HMM; each match
             state carries the region label of its alignment column
             (majority of the training precursors); residues take the label
             of the state they align to (inserts: the preceding state).
             Unaligned ends next to a partial alignment are extended
             ungapped over the remaining states: within the region of the
             last aligned state only (extension "within", so a mature
             peptide alone does not grow a propeptide), or across regions
             ("cross"); past the model ends, up to `overhang` residues take
             the terminal label.
  envelope   the best domain envelope of each region-profile HMM, put in
             SIG < PRO < MAT < POST order, overlaps split in the middle.
  cterm      optional C-terminal processing: a mature peptide that ends the
             sequence (nothing predicted after it) loses a trailing
             G[KR]{0,3}, the glycine that donates the C-terminal amide and
             the basic residues removed before amidation.
  sites      optional refinement of the three cleavage boundaries (signal
             peptide end, mature start, mature end) with lightmotif
             position weight matrices trained on exact UniProt boundaries:
             within +/- radius residues, the boundary maximising
             log-odds - penalty x |shift| is kept.

Database and features (used by conodictor; conodb validates them against a
release)
  Database   a conodb package directory (profiles.hmm/.prf/.pwm.json,
             precursors.hmm + labels, sites.json, framework.json,
             model.json, settings.json, manifest.json)
  compute    every feature of the model for every sequence, exactly as in
             training: HMM best-domain bits (pyhmmer, fixed Z, max mode),
             PSSM raw scores (pfsearch -f -k -a -r), best lightmotif window
             score, and the framework of the predicted mature peptide; then
             quantised to 4 significant digits like the values the model was
             fitted on (features.tsv).

This file is copied verbatim into conodictor (CORE_VERSION must match the
package's core_version; packages without one were built with version 1).
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

CORE_VERSION = 1

REGIONS = ("SIG", "PRO", "MAT", "POST")
AA = "ACDEFGHIKLMNPQRSTVWY"
SITE_TYPES = ("signal", "mature_start", "mature_end")
NONE = "NONE"

Segments = dict[str, tuple[int, int]]


def clean(seq: str) -> str:
    """Upper-case protein sequence; non-letters removed, '*' (stop) removed."""
    return re.sub(r"[^A-Z]", "", seq.upper())


def cysteine_framework(seq: str) -> str:
    """Cysteine pattern, e.g. 'GCCSDPRCAWRC' -> 'CC-C-C' ('' if no Cys)."""
    core = seq.upper().strip("ABDEFGHIJKLMNOPQRSTUVWXYZ")
    return re.sub(r"[^C]+", "-", core)


# ----------------------------------------------------------- HMM scoring ----

def _s(x) -> str:
    return x if isinstance(x, str) else x.decode()


@dataclass
class DomainHit:
    score: float
    start: int
    end: int
    hmm_from: int = 0
    hmm_to: int = 0
    hmm_aln: str = ""
    seq_aln: str = ""


def best_domains(hmms: list, seqs: Iterable[tuple[str, str]], *, min_bits: float = -1000.0, max_mode: bool = True,
                 cpus: int = 0, alignments: bool = False) -> dict[tuple[str, str], DomainHit]:
    """Best domain of every (sequence, HMM) pair HMMER reports: {(seq_id, hmm_name): DomainHit}.

    Fixed Z (E-values unused), optional max mode (no acceleration filters,
    so short profiles always score), reporting floor min_bits.
    """
    import os

    import pyhmmer
    from pyhmmer import easel

    seqs = list(seqs)
    if not seqs or not hmms:
        return {}
    alpha = easel.Alphabet.amino()
    block = easel.TextSequenceBlock(
        [easel.TextSequence(name=i.encode(), sequence=re.sub(r"[^A-Z]", "X", s.upper())) for i, s in seqs]
    ).digitize(alpha)
    kw = dict(Z=1, domZ=1, T=min_bits, domT=min_bits, cpus=cpus or os.cpu_count() or 1)
    if max_mode:
        kw |= dict(F1=1.0, F2=1.0, F3=1.0, bias_filter=False)
    out: dict[tuple[str, str], DomainHit] = {}
    for hits in pyhmmer.hmmsearch(hmms, block, **kw):
        q = _s(hits.query.name)
        for h in hits:
            doms = [d for d in h.domains if d.score >= min_bits]
            if not doms:
                continue
            d = max(doms, key=lambda d: (d.score, -d.env_from))
            hit = DomainHit(round(float(d.score), 2), d.env_from, d.env_to)
            if alignments:
                a = d.alignment
                hit.hmm_from, hit.hmm_to = a.hmm_from, a.hmm_to
                hit.start, hit.end = a.target_from, a.target_to
                hit.hmm_aln, hit.seq_aln = a.hmm_sequence, a.target_sequence
            out[(_s(h.name), q)] = hit
    return out


# ------------------------------------------------------ region prediction ----

def transfer_labels(seqlen: int, hit: DomainHit, labels: list[str], overhang: int = 5,
                    extension: str = "within") -> list[str]:
    """Per-residue region labels ('-' = none) from an alignment to a labelled full-precursor HMM.

    Unaligned ends are extended ungapped over the remaining model states,
    staying in the region of the last aligned state (extension "within") or
    crossing into other regions ("cross"). Where the extension stops (model
    end, or with "within" a region boundary), up to `overhang` further
    residues keep the label reached: precursor ends and peptide termini are
    often longer than the match states, which keep only columns occupied in
    most training sequences. In "cross" mode the extension only stops at
    the model ends, so the overhang applies there only."""
    if extension not in ("within", "cross"):
        raise ValueError(f"unknown extension {extension!r}")
    cross = extension == "cross"
    M = len(labels)
    out = ["-"] * seqlen
    k, t = hit.hmm_from, hit.start          # 1-based model state, sequence position
    for h, q in zip(hit.hmm_aln, hit.seq_aln, strict=True):
        if h == ".":                         # insert: residue between states
            if 1 <= t <= seqlen:
                out[t - 1] = labels[max(0, k - 2)] if k > 1 else labels[0]
            t += 1
        elif q == "-":                       # deletion: state without residue
            k += 1
        else:
            if 1 <= t <= seqlen and 1 <= k <= M:
                out[t - 1] = labels[k - 1]
            t += 1
            k += 1
    # N-terminal side
    i, s = hit.start - 1, hit.hmm_from - 1
    first = labels[min(M, hit.hmm_from) - 1]
    while i >= 1 and s >= 1 and (cross or labels[s - 1] == first):
        out[i - 1] = labels[s - 1]
        i, s = i - 1, s - 1
    tail = labels[0] if s == 0 else first    # model start reached, or (within) region boundary
    for j in range(i, max(0, i - overhang), -1):
        out[j - 1] = tail
    # C-terminal side
    i, s = hit.end + 1, hit.hmm_to + 1
    last = labels[min(M, hit.hmm_to) - 1]
    while i <= seqlen and s <= M and (cross or labels[s - 1] == last):
        out[i - 1] = labels[s - 1]
        i, s = i + 1, s + 1
    tail = labels[-1] if s == M + 1 else last
    for j in range(i, min(seqlen, i + overhang - 1) + 1):
        out[j - 1] = tail
    return out


def _runs(labels: list[str]) -> dict[str, tuple[int, int]]:
    """Longest contiguous run of each region label."""
    best: dict[str, tuple[int, int]] = {}
    i = 0
    while i < len(labels):
        j = i
        while j + 1 < len(labels) and labels[j + 1] == labels[i]:
            j += 1
        r = labels[i]
        if r in REGIONS and (r not in best or (j - i) > best[r][1] - best[r][0]):
            best[r] = (i + 1, j + 1)
        i = j + 1
    return best


def order_segments(seg: Segments, min_run: int) -> Segments:
    """Keep regions in SIG < PRO < MAT < POST order (the mature peptide has
    priority), split overlaps in the middle, drop runs shorter than min_run."""
    seg = {r: v for r, v in seg.items() if v[1] - v[0] + 1 >= min_run}
    keep: Segments = {}
    if "MAT" in seg:
        keep["MAT"] = seg["MAT"]
        ms, me = seg["MAT"]
        for r in ("SIG", "PRO"):
            if r in seg and seg[r][0] < ms:
                keep[r] = (seg[r][0], min(seg[r][1], ms - 1))
        if "POST" in seg and seg["POST"][1] > me:
            keep["POST"] = (max(seg["POST"][0], me + 1), seg["POST"][1])
    else:
        keep = dict(seg)
    order = [r for r in REGIONS if r in keep]
    for a, b in zip(order, order[1:], strict=False):
        if keep[a][1] >= keep[b][0]:
            mid = (keep[a][1] + keep[b][0]) // 2
            keep[a] = (keep[a][0], min(keep[a][1], mid))
            keep[b] = (max(keep[b][0], mid + 1), keep[b][1])
    return {r: v for r, v in keep.items() if v[1] - v[0] + 1 >= min_run}


def segments_from_labels(labels: list[str], min_run: int) -> Segments:
    return order_segments(_runs(labels), min_run)


def envelope_segments(region_hits: dict[str, tuple[int, int]], min_run: int) -> Segments:
    return order_segments(dict(region_hits), min_run)


# ----------------------------------------------------------------- sites ----

def site_boundaries(seg: Segments, seqlen: int) -> dict[str, int]:
    """Boundary positions b (last residue before the cut), where both sides exist."""
    out = {}
    if "SIG" in seg and seg["SIG"][1] < seqlen:
        out["signal"] = seg["SIG"][1]
    if "MAT" in seg and seg["MAT"][0] > 1:
        out["mature_start"] = seg["MAT"][0] - 1
    if "MAT" in seg and seg["MAT"][1] < seqlen:
        out["mature_end"] = seg["MAT"][1]
    return out


def site_window(seq: str, b: int, left: int, right: int) -> str | None:
    """Residues b-left+1 .. b+right (1-based), or None if out of range / non-standard."""
    if b - left < 0 or b + right > len(seq):
        return None
    w = seq[b - left:b + right].upper()
    return w if all(c in AA for c in w) else None


@dataclass
class SiteModels:
    left: int
    right: int
    matrices: dict  # site type -> lightmotif ScoringMatrix

    @classmethod
    def from_windows(cls, windows: dict[str, list[str]], left: int, right: int, pseudocount: float) -> SiteModels:
        import lightmotif
        mats = {}
        for t, ws in windows.items():
            ws = [w for w in ws if len(w) == left + right]
            if len(ws) >= 5:
                mats[t] = lightmotif.create(ws, protein=True).counts.normalize(pseudocount).log_odds()
        return cls(left, right, mats)

    def scores(self, seq: str, site: str) -> list[float] | None:
        """Log-odds of the window starting at each 0-based position."""
        import lightmotif
        m = self.matrices.get(site)
        L = self.left + self.right
        if m is None or len(seq) < L:
            return None
        s = re.sub(f"[^{AA}]", "A", seq.upper())
        striped = m.calculate(lightmotif.stripe(s, protein=True))
        return [float(striped[i]) for i in range(len(s) - L + 1)]


def refine_segments(seq: str, seg: Segments, sites: SiteModels, radius: int, penalty: float,
                    min_run: int, site_types: tuple[str, ...] | list[str] = SITE_TYPES) -> Segments:
    """Move each chosen cleavage boundary to the best-scoring site within +/- radius.
    A boundary is refined only where residues exist on both sides of it."""
    seg = dict(seg)
    for site in [t for t in SITE_TYPES if t in site_types]:
        b = site_boundaries(seg, len(seq)).get(site)
        if b is None:
            continue
        sc = sites.scores(seq, site)
        if sc is None:
            continue
        best, best_val = b, None
        for nb in range(b - radius, b + radius + 1):
            i = nb - sites.left
            if 0 <= i < len(sc):
                v = sc[i] - penalty * abs(nb - b)
                if best_val is None or v > best_val or (v == best_val and abs(nb - b) < abs(best - b)):
                    best, best_val = nb, v
        if best == b:
            continue
        if site == "signal":
            seg["SIG"] = (seg["SIG"][0], best)
            nxt = "PRO" if "PRO" in seg else ("MAT" if "MAT" in seg else None)
            if nxt and seg[nxt][0] <= b + 1 + radius:
                seg[nxt] = (best + 1, seg[nxt][1])
        elif site == "mature_start":
            seg["MAT"] = (best + 1, seg["MAT"][1])
            prv = "PRO" if "PRO" in seg else ("SIG" if "SIG" in seg else None)
            if prv and seg[prv][1] >= b - radius:
                seg[prv] = (seg[prv][0], best)
        else:
            seg["MAT"] = (seg["MAT"][0], best)
            if "POST" in seg and seg["POST"][0] <= b + 1 + radius:
                seg["POST"] = (best + 1, seg["POST"][1])
    return order_segments(seg, min_run)


_CTERM = re.compile(r"G[KR]{0,3}$")


def trim_cterm(seq: str, seg: Segments, min_run: int) -> tuple[Segments, bool]:
    """Remove a trailing G[KR]{0,3} from a mature peptide that nothing follows.
    Returns the segments and whether a C-terminal amidation is predicted."""
    if "MAT" not in seg or seg.get("POST"):
        return seg, False
    s, e = seg["MAT"]
    m = _CTERM.search(seq[s - 1:e].upper())
    if not m or (e - s + 1) - len(m.group()) < min_run:
        return seg, False
    out = dict(seg)
    out["MAT"] = (s, e - len(m.group()))
    return out, True


def predict_segments(seq: str, *, method: str, full_hit: DomainHit | None = None, labels: list[str] | None = None,
                     region_hits: dict[str, tuple[int, int]] | None = None, sites: SiteModels | None = None,
                     radius: int = 6, penalty: float = 0.5, min_run: int = 3, overhang: int = 5,
                     extension: str = "within", site_types: tuple[str, ...] | list[str] | None = None,
                     cterm: bool = False) -> Segments:
    """Regions of one sequence for one candidate family. `site_types`: which
    boundaries the site matrices refine (None = all three); `cterm`: apply
    the C-terminal G[KR]* processing rule (`trim_cterm`)."""
    if method == "transfer":
        if full_hit is None or not labels:
            return {}
        seg = segments_from_labels(transfer_labels(len(seq), full_hit, labels, overhang, extension), min_run)
    elif method == "envelope":
        seg = envelope_segments(region_hits or {}, min_run)
    else:
        raise ValueError(f"unknown region method {method!r}")
    types = SITE_TYPES if site_types is None else tuple(site_types)
    if sites is not None and seg and types:
        seg = refine_segments(seq, seg, sites, radius, penalty, min_run, types)
    if cterm and seg:
        seg = trim_cterm(seq, seg, min_run)[0]
    return seg


# ------------------------------------------------------------- framework ----

def mature_framework(seq: str, seg: Segments, short_peptide_length: int) -> str | None:
    """Cysteine framework of the predicted mature peptide ('-' = no cysteine);
    a short sequence without a prediction is taken as a mature peptide."""
    if "MAT" in seg:
        s, e = seg["MAT"]
        return cysteine_framework(seq[s - 1:e]) or "-"
    if len(seq) <= short_peptide_length:
        return cysteine_framework(seq) or "-"
    return None


def framework_logprob(fw: str | None, counts: Counter, vocab: int, pseudocount: float) -> float:
    """log P(framework | family), smoothed; None (no mature peptide) counts as unseen."""
    n = sum(counts.values())
    c = counts.get(fw, 0) if fw is not None else 0
    return math.log((c + pseudocount) / (n + pseudocount * vocab))


# ---------------------------------------------------------- model apply ----

def apply_model(model: dict, raw: np.ndarray) -> dict:
    """Class probabilities for one sequence.

    model: a model JSON written by conodb (families, features, center,
    scale, alpha, B, hmm_floor, fill, available). raw: K x P raw scores in
    model['features'] order, nan where the sequence got no score.
    Returns probabilities over families + NONE, P(NONE), and the best family
    with its probability conditional on the sequence being a conopeptide.
    """
    fams = model["families"]
    feats = model["features"]
    avail = np.array(model["available"], dtype=bool)
    center, scale = np.array(model["center"]), np.array(model["scale"])
    floor = model["hmm_floor"]
    fill = np.array([np.nan if v is None else v for v in model["fill"]])
    x = np.array(raw, dtype=float).copy()
    for p, f in enumerate(feats):
        col = x[:, p]
        if f.endswith(":hmm"):
            col = np.where(np.isnan(col), floor, np.maximum(col, floor))
        else:
            col = np.where(np.isnan(col), fill[p], col)
        x[:, p] = col
    z = np.where(avail, (x - center) / scale, 0.0)
    eta = np.array(model["alpha"]) + np.sum(z * np.array(model["B"]), axis=1)
    full = np.append(eta, 0.0)
    full = np.exp(full - full.max())
    prob = full / full.sum()
    k = int(np.argmax(prob[:-1]))
    fam_prob = prob[:-1] / max(prob[:-1].sum(), 1e-300)
    return {"probabilities": dict(zip(fams + [NONE], prob.round(6).tolist(), strict=True)),
            "prediction": (fams + [NONE])[int(np.argmax(prob))], "p_none": float(prob[-1]),
            "best_family": fams[k], "best_family_conditional": float(fam_prob[k])}


# -------------------------------------------------------------- features ----

def quantise(v) -> float:
    """A feature value as the model saw it in training: features.tsv stores
    floats with 4 significant digits and integers (PSSM raw scores) exactly."""
    if v is None:
        return float("nan")
    if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
        return float(v)
    v = float(v)
    return v if math.isnan(v) else float(f"{v:.4g}")


# xPSA header lines, with or without "/start-end", e.g.
#   >q0 motif=CONOPEP_A_MAT|CONOPEP_A_MAT raw_score=911
#   >q0/12-40 level=NA raw_score=911 seq_end=-5
_PSA = re.compile(r"^>(\S+?)(?:/(-?\d+)-(-?\d+))?\s.*?raw_score=\s*(-?\d+)")


class PfsearchError(RuntimeError):
    pass


def split_profiles(text: str) -> dict[str, str]:
    """PROSITE profile blocks of a .prf file by profile ID."""
    out = {}
    for blk in text.split("//\n"):
        if not blk.strip():
            continue
        m = re.search(r"^ID   (\S+);", blk, re.M)
        if m:
            out[m.group(1)] = blk + "//\n"
    return out


def pssm_best_scores(profile: str, seqs: list[tuple[str, str]], pfsearch: str,
                     options: Iterable[str] = ()) -> dict[str, int]:
    """Raw score of the best alignment of every sequence to one profile
    (text of a single PROSITE profile): pfsearch -f -k -a -r."""
    if not seqs:
        return {}
    with tempfile.TemporaryDirectory(prefix="conocore-ps-") as tmp:
        prf, fa = Path(tmp) / "p.prf", Path(tmp) / "s.fa"
        prf.write_text(profile)
        alias = {f"q{i}": sid for i, (sid, _) in enumerate(seqs)}      # pftools truncates long names
        fa.write_text("".join(f">q{i}\n{re.sub(r'[^A-Z]', 'X', s.upper())}\n" for i, (_, s) in enumerate(seqs)))
        cmd = [pfsearch, "-f", "-k", "-a", "-r", *options, str(prf), str(fa)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise PfsearchError(f"{' '.join(cmd)} failed (exit {res.returncode}): {res.stderr.strip()[:500]}")
    best: dict[str, int] = {}
    for line in res.stdout.splitlines():
        m = _PSA.match(line)
        if m and m.group(1) in alias:
            sid, sc = alias[m.group(1)], int(m.group(4))
            if sid not in best or sc > best[sid]:
                best[sid] = sc
    if len(best) != len(seqs):
        raise PfsearchError(f"pfsearch reported {len(best)} of {len(seqs)} sequences; output not understood")
    return best


def pwm_matrix(windows: list[str], pseudocount: float):
    import lightmotif
    return lightmotif.create(windows, protein=True).counts.normalize(pseudocount).log_odds()


def pwm_best(matrix, length: int, seq: str) -> float | None:
    """Best ungapped window log-odds (rounded to 3 decimals); None if the sequence is shorter."""
    import lightmotif
    q = re.sub(f"[^{AA}]", "A", seq.upper())
    if len(q) < length:
        return None
    sc = matrix.calculate(lightmotif.stripe(q, protein=True))
    return round(max(float(sc[i]) for i in range(len(q) - length + 1)), 3)


def _name(fam: str, region: str) -> str:
    return f"CONOPEP_{fam}_{region}"


@dataclass
class Database:
    """A conodb package, loaded for prediction."""
    path: Path
    model: dict
    settings: dict
    manifest: dict
    hmms: list
    profiles: dict[str, str]
    pwm_windows: dict[str, list[str]]
    full_hmms: list
    full_labels: dict[str, list[str]]
    sites: SiteModels | None
    framework: dict
    _pwm: dict = field(default_factory=dict, repr=False)

    @property
    def name(self) -> str:
        return self.manifest.get("database", self.path.name)

    @property
    def families(self) -> list[str]:
        return list(self.model["families"])

    @property
    def features(self) -> list[str]:
        return list(self.model["features"])

    def methods(self) -> set[str]:
        return {f.split(":")[1] for f in self.features}

    @classmethod
    def load(cls, path: str | Path) -> Database:
        from pyhmmer import plan7

        path = Path(path)
        need = ("model.json", "settings.json", "profiles.hmm")
        missing = [f for f in need if not (path / f).exists()]
        if missing:
            raise FileNotFoundError(f"{path} is not a conodb package (missing {', '.join(missing)})")

        def js(name, default=None):
            p = path / name
            return json.loads(p.read_text()) if p.exists() else default

        manifest = js("manifest.json", {})
        cv = int(manifest.get("core_version", 1))
        if cv != CORE_VERSION:
            raise ValueError(f"{path.name} needs core version {cv}, this program has {CORE_VERSION}")
        settings = js("settings.json")
        with plan7.HMMFile(path / "profiles.hmm") as fh:
            hmms = list(fh)
        full_hmms = []
        if (path / "precursors.hmm").exists() and (path / "precursors.hmm").stat().st_size:
            with plan7.HMMFile(path / "precursors.hmm") as fh:
                full_hmms = list(fh)
        an = settings.get("anatomy", {})
        sd = js("sites.json")
        sites = SiteModels.from_windows(sd["windows"], sd["left"], sd["right"], an.get("pseudocount", 0.5)) \
            if sd else None
        prf = (path / "profiles.prf").read_text() if (path / "profiles.prf").exists() else ""
        return cls(path=path, model=js("model.json"), settings=settings, manifest=manifest, hmms=hmms,
                   profiles=split_profiles(prf), pwm_windows=js("profiles.pwm.json", {}), full_hmms=full_hmms,
                   full_labels=js("precursors.labels.json", {}), sites=sites,
                   framework=js("framework.json", {"families": {}, "vocab_size": 1}))

    def region_kwargs(self) -> dict:
        an = self.settings.get("anatomy", {})
        s = an.get("sites", True)
        types = SITE_TYPES if s is True else () if s is False else tuple(t for t in SITE_TYPES if t in s)
        return {"method": an.get("method", "transfer"), "radius": an.get("search_radius", 6),
                "penalty": an.get("distance_penalty", 0.5), "min_run": an.get("min_run", 3),
                "overhang": an.get("terminal_overhang", 5), "extension": an.get("extension", "within"),
                "site_types": types}

    def pwm(self, name: str):
        if name not in self._pwm:
            w = self.pwm_windows.get(name) or []
            pc = self.settings.get("anatomy", {}).get("pseudocount", 0.5)
            self._pwm[name] = (pwm_matrix(w, pc), len(w[0])) if w else None
        return self._pwm[name]


def hmmer_candidates(db: Database, seqs: list[tuple[str, str]], *, f3: float = 1e-3, cpus: int = 0) -> set[str]:
    """Sequences that pass HMMER's acceleration filters (MSV F1 = 0.02, Viterbi
    F2 = 1e-3, Forward F3 = `f3`, per-comparison P-values) for at least one
    region or precursor HMM. A fast, optional prefilter: sequences outside
    this set are very unlikely to be conopeptides and are not scored."""
    import pyhmmer
    from pyhmmer import easel

    if not seqs:
        return set()
    block = easel.TextSequenceBlock(
        [easel.TextSequence(name=i.encode(), sequence=re.sub(r"[^A-Z]", "X", s.upper())) for i, s in seqs]
    ).digitize(easel.Alphabet.amino())
    out: set[str] = set()
    for hmms in (db.hmms, db.full_hmms):
        if not hmms:
            continue
        for hits in pyhmmer.hmmsearch(hmms, block, Z=1, domZ=1, T=-1000, domT=-1000, F3=f3,
                                      cpus=cpus or os.cpu_count() or 1):
            out.update(_s(h.name) for h in hits)
    return out


@dataclass
class SequenceResult:
    raw: np.ndarray                      # K x P quantised features, model order
    hmm_max: float                       # best HMM bits over the model's HMM features
    segments: dict[str, Segments]        # per family: predicted regions
    amidated: dict[str, bool]            # per family: C-terminal G[KR]* trimmed
    frameworks: dict[str, str | None]    # per family: framework of the predicted mature peptide


def compute_features(db: Database, seqs: list[tuple[str, str]], *, cpus: int = 0, pfsearch: str | None = None,
                     only: set[str] | None = None) -> dict[str, SequenceResult]:
    """Features of every sequence for db.model, plus the region predictions behind
    the framework feature. With `only` (e.g. from hmmer_candidates), other
    sequences are skipped entirely and absent from the result."""
    sc = db.settings.get("score", {})
    an = db.settings.get("anatomy", {})
    min_bits, max_mode = sc.get("hmm_min_bits", -1000.0), sc.get("max_mode", True)
    fams, feats = db.families, db.features
    K, P = len(fams), len(feats)
    fidx = {f: j for j, f in enumerate(feats)}
    used = [(f.split(":")[0], f.split(":")[1]) for f in feats]          # (region, method)
    hmm_regions = [rg for rg, m in used if m == "hmm"]
    seqs = [(i, clean(s)) for i, s in seqs if only is None or i in only]
    raw = {i: np.full((K, P), np.nan) for i, _ in seqs}

    hits = best_domains(db.hmms, seqs, min_bits=min_bits, max_mode=max_mode, cpus=cpus)
    hmm_max = {i: -math.inf for i, _ in seqs}
    for k, fam in enumerate(fams):
        for rg in hmm_regions:
            name, j = _name(fam, rg), fidx[f"{rg}:hmm"]
            for i, _ in seqs:
                h = hits.get((i, name))
                if h is not None:
                    raw[i][k, j] = quantise(h.score)
                    hmm_max[i] = max(hmm_max[i], h.score)
    todo = seqs

    if any(m == "pssm" for _, m in used) and todo:
        if not pfsearch:
            raise PfsearchError("the model uses PSSM features: pfsearch (pftools 2.3) is required")
        opts = sc.get("pfsearch_options", [])
        jobs = [(k, rg, db.profiles[_name(fam, rg)]) for k, fam in enumerate(fams) for rg, m in used
                if m == "pssm" and _name(fam, rg) in db.profiles]
        with ThreadPoolExecutor(max_workers=cpus or os.cpu_count() or 1) as pool:
            results = list(pool.map(lambda jb: pssm_best_scores(jb[2], todo, pfsearch, opts), jobs))
        for (k, rg, _), best in zip(jobs, results, strict=True):
            j = fidx[f"{rg}:pssm"]
            for i, v in best.items():
                raw[i][k, j] = quantise(int(v))

    for k, fam in enumerate(fams):
        for rg, m in used:
            if m != "pwm":
                continue
            mat = db.pwm(_name(fam, rg))
            if mat is None:
                continue
            j = fidx[f"{rg}:pwm"]
            for i, s in todo:
                v = pwm_best(mat[0], mat[1], s)
                if v is not None:
                    raw[i][k, j] = quantise(v)

    segs: dict[str, dict] = {i: {} for i, _ in seqs}
    amid: dict[str, dict] = {i: {} for i, _ in seqs}
    fws: dict[str, dict] = {i: {} for i, _ in seqs}
    if todo and db.full_labels:
        kw = db.region_kwargs()
        cterm = bool(an.get("cterm_trim", False))
        sites = db.sites if kw["site_types"] else None
        fhits = best_domains(db.full_hmms, todo, min_bits=min_bits, max_mode=max_mode, cpus=cpus,
                             alignments=True)
        fw_j = fidx.get("FW:framework")
        fwc = db.framework.get("families", {})
        V = int(db.framework.get("vocab_size", 1))
        pc = sc.get("framework_pseudocount", 0.5)
        short = an.get("short_peptide_length", 60)
        for k, fam in enumerate(fams):
            name = _name(fam, "FULL")
            for i, s in todo:
                seg = predict_segments(s, full_hit=fhits.get((i, name)), labels=db.full_labels.get(name),
                                       sites=sites, **kw)
                a = False
                if cterm and seg:
                    seg, a = trim_cterm(s, seg, kw["min_run"])
                fw = mature_framework(s, seg, short)
                segs[i][fam], amid[i][fam], fws[i][fam] = seg, a, fw
                if fw_j is not None:
                    raw[i][k, fw_j] = quantise(round(framework_logprob(fw, Counter(fwc.get(fam, {})), V, pc), 4))
    return {i: SequenceResult(raw[i], hmm_max[i], segs[i], amid[i], fws[i]) for i, _ in seqs}


def find_pfsearch(name: str = "pfsearch") -> str | None:
    return shutil.which(name)
