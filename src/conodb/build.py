"""Build: region profiles (HMM with pyhmmer, PSSM with pftools) per superfamily.

Profile sets:

  final     one profile per (superfamily, region) from all training entries
  holdout   the same family's profiles rebuilt without one cross-validation
            unit: each fold (tier 1) or each entry group (tier 2, LOGO)
  ship      the same from every labelled entry (train + test), for the
            database conodictor uses once evaluation is done

final and holdout use 'train' entries only (never test, temporal_test,
excluded or reference entries).

Per (superfamily, set, unit), regions SIG/PRO/MAT/POST give:
  .afa   FAMSA alignment     .hmm  pyhmmer HMM     .prf  pftools PSSM
  .pwm.json  ungapped lightmotif windows (columns occupied >= min_occupancy;
             gaps there filled with the column consensus)
and region FULL (complete precursors) gives a full-precursor HMM with hand
architecture: match states = columns occupied >= anatomy.min_occupancy, each
labelled with the majority UniProt region of its residues (.labels.json).
Per set, cleavage-site windows from exact UniProt boundaries (.sites.json).

Out-of-fold scoring then needs, for an entry, its own family's holdout
profiles for its unit plus the final profiles of every other family (which
never contain it). That gives leakage-free features for calibrating cutoffs
and fitting the multinomial model.

Every profile has a cache key (member sequences + build settings + tool
versions); a rebuild reuses any profile whose key and file are unchanged,
so a database update only rebuilds families whose data changed.

Output names follow ConoDictor: CONOPEP_<superfamily>_<REGION>.
"""

from __future__ import annotations

import glob
import io
import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import pyfamsa
import pyhmmer
from pyhmmer import easel, plan7

from .config import BuildConfig, Config
from .io import atomic_open, read_tsv, sha256_file, sha256_obj

log = logging.getLogger(__name__)

# Bump when the content of built files changes for identical inputs, so
# cached profiles from an older conodb are rebuilt.
BUILD_FORMAT = 4

PROFILE_FIELDS = [
    "set",
    "superfamily",
    "region",
    "unit",
    "name",
    "status",
    "n_sequences",
    "n_unique",
    "n_entries",
    "n_groups",
    "match_length",
    "hmm",
    "pssm",
    "msa",
    "labels",
    "pwm",
    "cache_key",
    "hmm_sha256",
    "pssm_sha256",
]
_LABEL_CODE = {"SIG": "S", "PRO": "P", "MAT": "M", "POST": "T", "INTER": "I"}
_CODE_LABEL = {v: k for k, v in _LABEL_CODE.items()} | {"-": "-"}


class BuildError(RuntimeError):
    pass


def safe(name: str) -> str:
    """Superfamily label -> token usable in profile names and paths."""
    return re.sub(r"[^A-Za-z0-9.+-]+", "-", name).strip("-") or "unnamed"


@dataclass
class Spec:
    set: str  # final | holdout | ship
    superfamily: str
    region: str  # SIG | PRO | MAT | POST | FULL
    unit: str  # "" for final/ship; "fold3" / "group-P01234" for holdout
    members: list[tuple[str, str]]  # (region_id or accession, sequence)
    entries: set[str] = field(default_factory=set)
    groups: set[str] = field(default_factory=set)
    labels: dict[str, str] = field(
        default_factory=dict
    )  # FULL: member -> per-residue region codes

    @property
    def name(self) -> str:
        return f"CONOPEP_{safe(self.superfamily)}_{self.region}"

    @property
    def stem(self) -> Path:
        unit = self.unit or "all"
        return Path("profiles") / self.set / safe(self.superfamily) / unit / self.region


#  pftools ----


@dataclass
class PfTools:
    pfw: str
    pfmake: str
    matrix: Path
    version: str

    @classmethod
    def locate(cls, bc: BuildConfig) -> PfTools | None:
        pfw, pfmake = shutil.which(bc.pfw), shutil.which(bc.pfmake)
        if not pfw or not pfmake:
            return None
        matrix = _find_matrix(bc.matrix, Path(pfmake))
        if matrix is None:
            raise BuildError(
                f"pftools found but the score matrix {bc.matrix!r} was not. Set build.matrix to its full path "
                "(it ships with pftools, e.g. .../share/pftools*/blosum45.cmp or data/Matrices/blosum45.cmp)."
            )
        run = subprocess.run(
            [pfmake, "-h"], capture_output=True, text=True, check=False
        )
        first = next(
            (ln.strip() for ln in (run.stdout + run.stderr).splitlines() if ln.strip()),
            "pfmake ?",
        )
        return cls(pfw=pfw, pfmake=pfmake, matrix=matrix, version=first)


def _find_matrix(matrix: str, pfmake: Path) -> Path | None:
    p = Path(matrix).expanduser()
    if p.is_file():
        return p.resolve()
    name = p.name
    roots = [
        os.environ.get("PFTOOLS_DATA", ""),
        os.environ.get("CONDA_PREFIX", "")
        and str(Path(os.environ["CONDA_PREFIX"]) / "share"),
        str(pfmake.resolve().parent.parent / "share"),
        "/usr/local/share",
        "/usr/share",
        "/opt",
    ]
    for root in filter(None, roots):
        hits = sorted(glob.glob(os.path.join(root, "**", name), recursive=True))
        if hits:
            return Path(hits[0])
    return None


def _run(cmd: list[str]) -> str:
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0 or not res.stdout.strip():
        raise BuildError(
            f"{' '.join(cmd)} failed (exit {res.returncode}): {res.stderr.strip()[:500]}"
        )
    return res.stdout


def build_pssm(
    spec: Spec, aligned: list[tuple[str, str]], tools: PfTools, bc: BuildConfig
) -> str:
    """Aligned FASTA -> pfw (sequence weights) -> pfmake -> PROSITE profile text."""
    with tempfile.TemporaryDirectory(prefix="conodb-pf-") as tmp:
        afa = Path(tmp) / "in.msa"
        # short identifiers: pftools truncates long ones, and only weights matter here
        afa.write_text("".join(f">s{i}\n{seq}\n" for i, (_, seq) in enumerate(aligned)))
        weighted = Path(tmp) / "weighted.msa"
        weighted.write_text(_run([tools.pfw, "-m", *bc.pfw_options, str(afa)]))
        prf = _run(
            [tools.pfmake, "-m", *bc.pfmake_options, str(weighted), str(tools.matrix)]
        )
    return _label_profile(prf, spec)


def _label_profile(prf: str, spec: Spec) -> str:
    """Give the profile ConoDictor's identifiers and a reproducible header."""
    out, seen_ac = [], False
    for line in prf.splitlines():
        if line.startswith("ID   "):
            out.append(f"ID   {spec.name}; MATRIX.")
            continue
        if line.startswith("AC   "):
            out.append(f"AC   {spec.name};")
            seen_ac = True
            continue
        if line.startswith("DT   "):
            out.append("DT   built by conodb.")
            continue
        if line.startswith("CC   /GENERATED_BY="):
            # pfmake records its command line, including a random temp dir:
            # replace it so identical inputs give identical files.
            out.append('CC   /GENERATED_BY="conodb: pfw -m | pfmake -m";')
            continue
        if line.startswith("DE   "):
            out.append(
                f"DE   ConoDictor {spec.superfamily} {spec.region} region"
                + (f" (holdout {spec.unit})" if spec.unit else "")
                + f"; {len(spec.members)} sequences."
            )
            continue
        out.append(line)
    if not seen_ac:
        i = next((k for k, ln in enumerate(out) if ln.startswith("ID   ")), 0)
        out.insert(i + 1, f"AC   {spec.name};")
    return "\n".join(out).rstrip("\n") + "\n"


# ------------------------------------------------------------------- hmm ----


def align(members: list[tuple[str, str]], threads: int) -> list[tuple[str, str]]:
    seqs = [pyfamsa.Sequence(i.encode(), s.encode()) for i, s in members]
    msa = pyfamsa.Aligner(guide_tree="upgma", threads=threads).align(seqs)
    out = [(s.id.decode(), s.sequence.decode()) for s in msa]
    order = {i: k for k, (i, _) in enumerate(members)}
    return sorted(out, key=lambda x: order[x[0]])


def build_hmm(spec: Spec, aligned: list[tuple[str, str]]) -> tuple[bytes, int]:
    alpha = easel.Alphabet.amino()
    tm = easel.TextMSA(
        name=spec.name.encode(),
        sequences=[easel.TextSequence(name=i.encode(), sequence=s) for i, s in aligned],
    )
    hmm, _, _ = plan7.Builder(alpha, seed=42).build_msa(
        tm.digitize(alpha), plan7.Background(alpha)
    )
    hmm.creation_time = None  # reproducible files: no build date
    hmm.command_line = None
    desc = f"ConoDictor {spec.superfamily} {spec.region} region" + (
        f" (holdout {spec.unit})" if spec.unit else ""
    )
    hmm.description = desc.encode() if isinstance(hmm.name, bytes) else desc
    buf = io.BytesIO()
    hmm.write(buf)
    return buf.getvalue(), hmm.M


# ------------------------------------------------------- full / pwm / sites --


def column_occupancy(aligned: list[tuple[str, str]]) -> list[float]:
    n = len(aligned)
    return [sum(s[c] != "-" for _, s in aligned) / n for c in range(len(aligned[0][1]))]


def build_full(
    spec: Spec, aligned: list[tuple[str, str]], min_occupancy: float
) -> tuple[bytes, int, list[str]]:
    """Hand-architecture HMM whose match states are the well-occupied columns,
    each labelled with the majority region of the residues in that column."""
    occ = column_occupancy(aligned)
    match = [c for c, o in enumerate(occ) if o >= min_occupancy]
    labels = []
    for c in match:
        votes = []
        for mid, row in aligned:
            if row[c] != "-":
                idx = len(row[:c].replace("-", ""))
                votes.append(_CODE_LABEL.get(spec.labels[mid][idx], "-"))
        cnt = defaultdict(int)
        for v in votes:
            cnt[v] += 1
        labels.append(sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0])
    alpha = easel.Alphabet.amino()
    tm = easel.TextMSA(
        name=spec.name.encode(),
        sequences=[easel.TextSequence(name=i.encode(), sequence=r) for i, r in aligned],
    )
    tm.reference = "".join("x" if c in set(match) else "." for c in range(len(occ)))
    hmm, _, _ = plan7.Builder(alpha, seed=42, architecture="hand").build_msa(
        tm.digitize(alpha), plan7.Background(alpha)
    )
    hmm.creation_time = None
    hmm.command_line = None
    desc = f"ConoDictor {spec.superfamily} full precursor" + (
        f" (holdout {spec.unit})" if spec.unit else ""
    )
    hmm.description = desc.encode() if isinstance(hmm.name, bytes) else desc
    if hmm.M != len(labels):
        raise BuildError(
            f"{spec.name}: {hmm.M} match states but {len(labels)} labelled columns"
        )
    buf = io.BytesIO()
    hmm.write(buf)
    return buf.getvalue(), hmm.M, labels


def pwm_windows(aligned: list[tuple[str, str]], min_occupancy: float) -> list[str]:
    """Ungapped windows over the well-occupied columns (gaps -> column consensus)."""
    AA = "ACDEFGHIKLMNPQRSTVWY"
    occ = column_occupancy(aligned)
    cols = [c for c, o in enumerate(occ) if o >= min_occupancy]
    cons = {}
    for c in cols:
        cnt = defaultdict(int)
        for _, s in aligned:
            if s[c] in AA:
                cnt[s[c]] += 1
        cons[c] = (
            sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] if cnt else "A"
        )
    return (
        ["".join(s[c] if s[c] in AA else cons[c] for c in cols) for _, s in aligned]
        if cols
        else []
    )


def site_training_windows(
    entries_ok: set[str],
    regions: list[dict],
    seqs: dict[str, str],
    left: int,
    right: int,
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Cleavage-site windows from exact UniProt boundaries of the given entries,
    and the accession each window comes from (so cross-validation can leave
    a unit's own windows out)."""
    from .core import site_boundaries, site_window

    by_acc = defaultdict(dict)
    exact = defaultdict(lambda: True)
    for r in regions:
        if (
            r["accession"] in entries_ok
            and r["region"] in ("SIG", "PRO", "MAT", "POST")
            and r["index"] == "1"
        ):
            by_acc[r["accession"]][r["region"]] = (int(r["start"]), int(r["end"]))
            exact[(r["accession"], r["region"])] = r["exact"] == "true"
    out, accs = defaultdict(list), defaultdict(list)
    for acc, seg in sorted(by_acc.items()):
        seq = seqs.get(acc, "")
        for site, b in site_boundaries(seg, len(seq)).items():
            reg = {"signal": "SIG", "mature_start": "MAT", "mature_end": "MAT"}[site]
            if not exact[(acc, reg)]:
                continue
            w = site_window(seq, b, left, right)
            if w:
                out[site].append(w)
                accs[site].append(acc)
    return dict(out), dict(accs)


# ------------------------------------------------------------------ plan ----


def plan_specs(
    cfg: Config, entries: list[dict], regions: list[dict], partition: list[dict]
) -> tuple[list[Spec], list[str]]:
    bc = cfg.build
    part = {r["accession"]: r for r in partition}
    classes = {
        a
        for a, r in part.items()
        if r["tier"] in ("1", "2") and r["scope"] == "conopeptide"
    }
    train = {a for a in classes if part[a]["role"] == "train"}
    labelled = {
        a for a in classes if part[a]["role"] in ("train", "test", "temporal_test")
    }
    by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    dropped = defaultdict(int)
    for r in regions:
        if r["accession"] not in labelled or r["region"] not in bc.regions:
            continue
        if int(r["length"]) < bc.min_length:
            dropped[f"{r['region']} shorter than {bc.min_length} aa"] += 1
            continue
        if r["exact"] != "true" and not bc.use_uncertain:
            dropped[f"{r['region']} with uncertain boundaries"] += 1
            continue
        by_key[(part[r["accession"]]["superfamily"], r["region"])].append(r)

    # full precursors with per-residue region codes
    full_rows: dict[str, list[dict]] = defaultdict(list)
    codes: dict[str, str] = {}
    if bc.full_models:
        ent = {e["accession"]: e for e in entries}
        regs = defaultdict(list)
        for r in regions:
            regs[r["accession"]].append(r)
        for acc in sorted(labelled):
            e = ent.get(acc)
            if not e or e["complete_precursor"] != "true":
                continue
            lab = ["-"] * len(e["sequence"])
            for r in regs[acc]:
                for i in range(int(r["start"]) - 1, int(r["end"])):
                    lab[i] = _LABEL_CODE.get(r["region"], "-")
            codes[acc] = "".join(lab)
            full_rows[part[acc]["superfamily"]].append(
                {"region_id": acc, "accession": acc, "sequence": e["sequence"]}
            )
        for sf, rows in full_rows.items():
            by_key[(sf, "FULL")] = rows

    specs: list[Spec] = []
    for (sf, region), rows in sorted(by_key.items()):
        tier = {part[r["accession"]]["tier"] for r in rows}.pop()

        def make(set_name, unit, keep, sf=sf, region=region, rows=rows):
            sel = [r for r in rows if keep(r["accession"])]
            return Spec(
                set_name,
                sf,
                region,
                unit,
                [(r["region_id"], r["sequence"].upper()) for r in sel],
                {r["accession"] for r in sel},
                {part[r["accession"]]["group"] for r in sel},
                {r["region_id"]: codes[r["accession"]] for r in sel}
                if region == "FULL"
                else {},
            )

        specs.append(make("final", "", lambda a: a in train))
        specs.append(make("ship", "", lambda a: a in labelled))
        if tier == "1":
            for f in range(1, cfg.partition.folds + 1):
                specs.append(
                    make(
                        "holdout",
                        f"fold{f}",
                        lambda a, f=f: a in train and part[a]["fold"] != str(f),
                    )
                )
        else:
            for g in sorted(
                {part[r["accession"]]["group"] for r in rows if r["accession"] in train}
            ):
                specs.append(
                    make(
                        "holdout",
                        f"group-{safe(g)}",
                        lambda a, g=g: a in train and part[a]["group"] != g,
                    )
                )
    specs = [s for s in specs if s.members]
    notes = [f"{k}: {v} regions not used" for k, v in sorted(dropped.items())]
    return specs, notes


def _dedupe(members: list[tuple[str, str]]) -> list[tuple[str, str]]:
    seen: dict[str, str] = {}
    for rid, seq in sorted(members):
        seen.setdefault(seq, rid)
    return sorted(((rid, seq) for seq, rid in seen.items()), key=lambda x: x[0])


# ----------------------------------------------------------------- build ----


def run_build(
    cfg: Config,
    rdir: Path,
    entries: list[dict],
    regions: list[dict],
    partition: list[dict],
) -> tuple[list[dict], dict]:
    import json

    bc, an = cfg.build, cfg.anatomy
    tools = None
    if bc.pssm in ("auto", "true"):
        tools = PfTools.locate(bc)
        if tools is None and bc.pssm == "true":
            raise BuildError(
                f"build.pssm = true but {bc.pfw!r}/{bc.pfmake!r} not found on PATH"
            )
        if tools is None:
            log.warning(
                "pftools (pfw, pfmake) not found: building HMMs only (build.pssm = auto)"
            )
    use_pwm = bc.pwm
    if use_pwm:
        try:
            import lightmotif  # noqa: F401
        except ImportError:
            log.warning("lightmotif not installed: no PWMs (pip install lightmotif)")
            use_pwm = False
    specs, notes = plan_specs(cfg, entries, regions, partition)

    # leakage guard: final/holdout only training entries; holdouts never contain their unit
    part = {r["accession"]: r for r in partition}
    for s in specs:
        if s.set in ("final", "holdout"):
            bad = [a for a in s.entries if part[a]["role"] != "train"]
            if bad:
                raise BuildError(
                    f"{s.name} {s.set} would include non-training entries: {bad[:5]}"
                )
        else:
            bad = [
                a
                for a in s.entries
                if part[a]["role"] not in ("train", "test", "temporal_test")
            ]
            if bad:
                raise BuildError(
                    f"{s.name} ship would include unlabelled or excluded entries: {bad[:5]}"
                )
        if s.unit.startswith("fold") and any(
            part[a]["fold"] == s.unit[4:] for a in s.entries
        ):
            raise BuildError(f"{s.name} {s.unit} contains its own fold")
        if s.unit.startswith("group-") and any(
            f"group-{safe(part[a]['group'])}" == s.unit for a in s.entries
        ):
            raise BuildError(f"{s.name} {s.unit} contains its own group")

    previous = {}
    if (rdir / "profiles.tsv").exists():
        previous = {
            (r["set"], r["superfamily"], r["region"], r["unit"]): r
            for r in read_tsv(rdir / "profiles.tsv")
        }
    settings = {
        "format": BUILD_FORMAT,
        "build": cfg.section("build"),
        "occupancy": an.min_occupancy,
        "pyhmmer": pyhmmer.__version__,
        "pyfamsa": pyfamsa.__version__,
        "pftools": tools and tools.version,
        "matrix": tools and sha256_file(tools.matrix),
        "pwm": use_pwm,
    }
    threads = bc.threads or os.cpu_count() or 1

    def one(spec: Spec) -> dict:
        uniq = _dedupe(spec.members)
        row = {
            "set": spec.set,
            "superfamily": spec.superfamily,
            "region": spec.region,
            "unit": spec.unit,
            "name": spec.name,
            "n_sequences": len(spec.members),
            "n_unique": len(uniq),
            "n_entries": len(spec.entries),
            "n_groups": len(spec.groups),
        }
        key = sha256_obj(
            {
                "members": uniq,
                "labels": {i: spec.labels.get(i) for i, _ in uniq},
                "settings": settings,
            }
        )
        row["cache_key"] = key
        if len(uniq) < bc.min_unique:
            return row | {"status": f"skipped: {len(uniq)} unique < {bc.min_unique}"}
        stem = spec.stem
        hmm_p, pssm_p, msa_p = (
            stem.with_suffix(".hmm"),
            stem.with_suffix(".prf"),
            stem.with_suffix(".afa"),
        )
        lab_p, pwm_p = Path(f"{stem}.labels.json"), Path(f"{stem}.pwm.json")
        old = previous.get((spec.set, spec.superfamily, spec.region, spec.unit))
        if (
            old
            and old["cache_key"] == key
            and old["status"] == "built"
            and _intact(rdir, old)
        ):
            return old | {"status": "built"}
        aligned = align(uniq, 1)
        with atomic_open(rdir / msa_p) as fh:
            fh.write("".join(f">{i}\n{s}\n" for i, s in aligned))
        if spec.region == "FULL":
            hmm_bytes, m, labels = build_full(spec, aligned, an.min_occupancy)
            with atomic_open(rdir / lab_p) as fh:
                json.dump({"name": spec.name, "labels": labels}, fh)
            row |= {"labels": str(lab_p)}
        else:
            hmm_bytes, m = build_hmm(spec, aligned)
            if use_pwm:
                with atomic_open(rdir / pwm_p) as fh:
                    json.dump(
                        {
                            "name": spec.name,
                            "windows": pwm_windows(aligned, an.min_occupancy),
                        },
                        fh,
                    )
                row |= {"pwm": str(pwm_p)}
        with atomic_open(rdir / hmm_p, "wb") as fh:
            fh.write(hmm_bytes)
        row |= {
            "status": "built",
            "match_length": m,
            "hmm": str(hmm_p),
            "msa": str(msa_p),
            "hmm_sha256": sha256_file(rdir / hmm_p),
        }
        if tools is not None and spec.region != "FULL":
            with atomic_open(rdir / pssm_p) as fh:
                fh.write(build_pssm(spec, aligned, tools, bc))
            row |= {"pssm": str(pssm_p), "pssm_sha256": sha256_file(rdir / pssm_p)}
        return row

    with ThreadPoolExecutor(max_workers=threads) as pool:
        rows = list(pool.map(one, specs))
    reused = sum(
        1
        for r, s in zip(rows, specs, strict=True)
        if r.get("hmm")
        and previous.get((s.set, s.superfamily, s.region, s.unit), {}).get("cache_key")
        == r["cache_key"]
    )

    # libraries: regions (classification), full precursors (+labels), pwms
    written: set[Path] = set()
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["status"] != "built":
            continue
        stem = (
            f"profiles/{r['set']}"
            if r["set"] != "holdout"
            else f"profiles/holdout/{safe(r['superfamily'])}/{r['unit']}"
        )
        groups[stem].append(r)
    for stem in ("profiles/final", "profiles/ship"):  # always present, even if empty
        groups.setdefault(stem, [])
    for stem, rs in sorted(groups.items()):
        written |= _concat(rdir, [r for r in rs if r["region"] != "FULL"], stem)
        written |= _concat_full(rdir, [r for r in rs if r["region"] == "FULL"], stem)
        written |= _concat_pwm(rdir, [r for r in rs if r.get("pwm")], stem)

    # cleavage-site windows per set
    seqs = {e["accession"]: e["sequence"] for e in entries}
    for set_name, roles in (
        ("final", ("train",)),
        ("ship", ("train", "test", "temporal_test")),
    ):
        ok = {
            a
            for a, p in part.items()
            if p["role"] in roles
            and p["scope"] == "conopeptide"
            and p["tier"] in ("1", "2")
        }
        w, wa = site_training_windows(ok, regions, seqs, an.site_left, an.site_right)
        path = rdir / f"profiles/{set_name}.sites.json"
        with atomic_open(path) as fh:
            json.dump(
                {
                    "left": an.site_left,
                    "right": an.site_right,
                    "windows": w,
                    "accessions": wa,
                },
                fh,
                sort_keys=True,
            )
        written.add(path)
    _prune(rdir, rows, written)

    rows.sort(
        key=lambda r: (
            ["final", "ship", "holdout"].index(r["set"]),
            r["superfamily"],
            r["unit"],
            r["region"],
        )
    )
    built = [r for r in rows if r["status"] == "built"]
    stats = {
        "profiles": len(built),
        "final": sum(r["set"] == "final" for r in built),
        "ship": sum(r["set"] == "ship" for r in built),
        "holdout": sum(r["set"] == "holdout" for r in built),
        "full_models": sum(r["region"] == "FULL" for r in built),
        "reused": reused,
        "skipped": len(rows) - len(built),
        "pssm": tools is not None,
        "pwm": use_pwm,
        "notes": notes,
    }
    return rows, stats


def _intact(rdir: Path, row: dict) -> bool:
    for col, sha in (("hmm", "hmm_sha256"), ("pssm", "pssm_sha256")):
        if row.get(col):
            p = rdir / row[col]
            if not p.exists() or sha256_file(p) != row[sha]:
                return False
    return all((rdir / row[c]).exists() for c in ("labels", "pwm", "msa") if row.get(c))


_REGION_ORDER = {"SIG": 0, "PRO": 1, "MAT": 2, "POST": 3}


def _concat(rdir: Path, rows: list[dict], stem: str) -> set[Path]:
    rows = sorted(
        rows, key=lambda r: (r["superfamily"], _REGION_ORDER.get(r["region"], 9))
    )
    out = set()
    p = rdir / f"{stem}.hmm"
    with atomic_open(p, "wb") as fh:
        for r in rows:
            fh.write((rdir / r["hmm"]).read_bytes())
    out.add(p)
    if rows and all(r.get("pssm") for r in rows):
        p = rdir / f"{stem}.prf"
        with atomic_open(p) as fh:
            for r in rows:
                fh.write((rdir / r["pssm"]).read_text())
        out.add(p)
    return out


def _concat_full(rdir: Path, rows: list[dict], stem: str) -> set[Path]:
    import json

    if not rows:
        return set()
    rows = sorted(rows, key=lambda r: r["superfamily"])
    hp, lp = rdir / f"{stem}.full.hmm", rdir / f"{stem}.full.labels.json"
    with atomic_open(hp, "wb") as fh:
        for r in rows:
            fh.write((rdir / r["hmm"]).read_bytes())
    labels = {}
    for r in rows:
        with open(rdir / r["labels"]) as fh:
            d = json.load(fh)
        labels[d["name"]] = d["labels"]
    with atomic_open(lp) as fh:
        json.dump(labels, fh, sort_keys=True)
    return {hp, lp}


def _concat_pwm(rdir: Path, rows: list[dict], stem: str) -> set[Path]:
    import json

    if not rows:
        return set()
    lib = {}
    for r in sorted(
        rows, key=lambda r: (r["superfamily"], _REGION_ORDER.get(r["region"], 9))
    ):
        with open(rdir / r["pwm"]) as fh:
            d = json.load(fh)
        lib[d["name"]] = d["windows"]
    p = rdir / f"{stem}.pwm.json"
    with atomic_open(p) as fh:
        json.dump(lib, fh, sort_keys=True)
    return {p}


def _prune(rdir: Path, rows: list[dict], written: set[Path]) -> None:
    """Remove profile files no longer produced (e.g. a family lost its data)."""
    keep = {
        rdir / r[c]
        for r in rows
        for c in ("hmm", "pssm", "msa", "labels", "pwm")
        if r.get(c)
    } | written
    root = rdir / "profiles"
    if not root.exists():
        return
    for p in root.rglob("*"):
        if (
            p.is_file()
            and p.suffix in (".hmm", ".prf", ".afa", ".json")
            and p not in keep
        ):
            p.unlink()
