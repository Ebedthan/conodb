"""The variant grids compared by the select stage (cross-validation on
training rows) and reported on the test set by evaluate and anatomy.

One definition, so the cross-validation table used to choose the defaults and
the test-set table reported afterwards list exactly the same variants.

Classifier variants: HMM floor {censored at 0 bits, -20 bits} x feature types
{HMM, HMM+PSSM, HMM+PWM, HMM+PSSM+PWM} x {without, with the framework
feature}, regions SIG/PRO/MAT; plus the settings of ConoDictor's previous
database (censored, HMM+PSSM, POST included), signal only, without signal,
and the current settings if they are not already in the grid.

Region variants: envelope (with / without all sites, with / without the
C-terminal G[KR]* rule), and label transfer with extension {within, cross} x
overhang {0, 5} x sites {none, signal, mature start, signal + mature start,
all three} x {without, with the C-terminal rule}; plus the current settings
if not in the grid.

`complexity` orders variants for the selection rule (the simplest variant
whose cross-validated performance is not clearly worse than the best wins):
classifier = features per family; regions = boundaries refined + 1 if an
overhang is used + 1 for the C-terminal rule.
"""

from __future__ import annotations

from .config import Config

GRID_REGIONS = ("SIG", "PRO", "MAT")
FEATURE_SETS = (("hmm",), ("hmm", "pssm"), ("hmm", "pwm"), ("hmm", "pssm", "pwm"))
LABEL = {"hmm": "HMM", "pssm": "PSSM", "pwm": "PWM", "framework": "FW"}
SITE_ABBR = {"signal": "S", "mature_start": "MS", "mature_end": "ME"}
SITE_SETS = (
    (),
    ("signal",),
    ("mature_start",),
    ("signal", "mature_start"),
    ("signal", "mature_start", "mature_end"),
)


def _floor_name(f: float) -> str:
    return "censored" if f == 0 else f"floor {f:g}"


def classifier_name(methods, regions, floor) -> str:
    name = f"{_floor_name(floor)} · " + "+".join(LABEL[m] for m in methods)
    if tuple(regions) != GRID_REGIONS:
        name += " · " + "/".join(regions)
    return name


def _complexity(methods, regions) -> int:
    return sum(len(regions) for m in methods if m != "framework") + (
        "framework" in methods
    )


def classifier_variants(cfg: Config, have: set[str]) -> list[dict]:
    """Variants whose feature types were all scored (`have`), deduplicated."""
    out, seen = [], set()

    def add(methods, regions, floor, note=""):
        methods = tuple(m for m in ("hmm", "pssm", "pwm", "framework") if m in methods)
        regions = tuple(r for r in ("SIG", "PRO", "MAT", "POST") if r in regions)
        key = (methods, regions, float(floor))
        if key in seen or not set(methods) <= have or not methods or not regions:
            for v in out:  # keep the note on the existing entry
                if (
                    (v["methods"], v["regions"], v["hmm_floor"]) == key
                    and note
                    and note not in v["note"]
                ):
                    v["note"] = (v["note"] + "; " + note).strip("; ")
            return
        seen.add(key)
        out.append(
            {
                "name": classifier_name(methods, regions, floor),
                "methods": methods,
                "regions": regions,
                "hmm_floor": float(floor),
                "complexity": _complexity(methods, regions),
                "note": note,
            }
        )

    for floor in (0.0, -20.0):
        for fs in FEATURE_SETS:
            for fw in (False, True):
                add(fs + (("framework",) if fw else ()), GRID_REGIONS, floor)
    add(
        ("hmm", "pssm"),
        ("SIG", "PRO", "MAT", "POST"),
        0.0,
        "previous ConoDictor database",
    )
    fc = cfg.fit
    cur = tuple(m for m in fc.methods if m in have)
    add(cur, ("SIG",), fc.hmm_floor, "signal only")
    add(cur, tuple(r for r in fc.regions if r != "SIG"), fc.hmm_floor, "without signal")
    add(cur, tuple(fc.regions), fc.hmm_floor, "current settings")
    return out


def region_name(method, extension, overhang, sites, cterm=False) -> str:
    base = "envelope" if method == "envelope" else f"transfer/{extension}/oh{overhang}"
    return (
        base
        + (" + sites " + "+".join(SITE_ABBR[s] for s in sites) if sites else "")
        + (" + Cterm" if cterm else "")
    )


def region_variants(cfg: Config) -> list[dict]:
    out, seen = [], set()

    def add(method, extension, overhang, sites, cterm=False, note=""):
        sites = tuple(s for s in ("signal", "mature_start", "mature_end") if s in sites)
        if method == "envelope":
            extension, overhang = "-", 0
        key = (method, extension, int(overhang), sites, bool(cterm))
        if key in seen:
            for v in out:
                if v["key"] == key and note and note not in v["note"]:
                    v["note"] = (v["note"] + "; " + note).strip("; ")
            return
        seen.add(key)
        out.append(
            {
                "name": region_name(method, extension, overhang, sites, cterm),
                "key": key,
                "method": method,
                "extension": extension,
                "overhang": int(overhang),
                "sites": sites,
                "cterm": bool(cterm),
                "complexity": len(sites) + (overhang > 0) + bool(cterm),
                "note": note,
            }
        )

    for ct in (False, True):
        add("envelope", "-", 0, (), ct)
        add("envelope", "-", 0, SITE_SETS[-1], ct)
        for ext in ("within", "cross"):
            for oh in (0, 5):
                for ss in SITE_SETS:
                    add("transfer", ext, oh, ss, ct)
    an = cfg.anatomy
    add(
        an.method,
        an.extension,
        an.terminal_overhang,
        an.site_types,
        an.cterm_trim,
        "current settings",
    )
    return out


def region_kwargs(cfg: Config, v: dict) -> dict:
    """core.predict_segments keyword arguments of a region variant."""
    an = cfg.anatomy
    return {
        "method": v["method"],
        "radius": an.search_radius,
        "penalty": an.distance_penalty,
        "min_run": an.min_run,
        "overhang": v["overhang"],
        "extension": v["extension"] if v["method"] == "transfer" else "within",
        "site_types": v["sites"],
        "cterm": v["cterm"],
    }
