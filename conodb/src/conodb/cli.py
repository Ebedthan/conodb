"""Command line: conodb init | run | check | status | stages."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import EXAMPLE_CONFIG, EXAMPLE_LABELS, load_config
from .io import read_json
from .build import BuildError
from .fit import FitError
from .score import ScoreError
from .pipeline import STAGES, STAGE_NAMES, CurationNeeded, latest_release, run
from .uniprot import FetchError, probe

# Exit codes (useful for cron / CI)
OK, ERROR, CURATION, UPDATE_AVAILABLE = 0, 1, 2, 10


def _cmd_init(args) -> int:
    d = Path(args.directory)
    d.mkdir(parents=True, exist_ok=True)
    for name, text in (("conodb.toml", EXAMPLE_CONFIG), ("labels.tsv", EXAMPLE_LABELS)):
        p = d / name
        if p.exists():
            print(f"exists, left unchanged: {p}")
        else:
            p.write_text(text, encoding="utf-8")
            print(f"created: {p}")
    print(f"\nNext: conodb run -c {d / 'conodb.toml'}")
    return OK


def _cmd_run(args) -> int:
    cfg = load_config(args.config)
    if args.input and not args.input.exists():
        logging.error("Input file not found: %s", args.input)
        return ERROR
    rdir = run(cfg, force=args.force, from_stage=args.from_stage, until_stage=args.until,
               local_input=args.input, release_name=args.release)
    print(rdir / "report.md")
    return OK


def _cmd_check(args) -> int:
    cfg = load_config(args.config)
    current = latest_release(cfg.output.root)
    info = probe(cfg.source)
    if current == info.release:
        print(f"Up to date: UniProt {info.release} already built.")
        return OK
    print(f"Update available: UniProt {info.release} (built: {current or 'none'}), {info.total} entries.")
    return UPDATE_AVAILABLE


def _cmd_status(args) -> int:
    cfg = load_config(args.config)
    root = cfg.output.root
    name = latest_release(root)
    if name is None:
        print(f"No release built yet in {root}.")
        return OK
    m = read_json(root / "releases" / name / "manifest.json")
    print(f"Latest release: {name} (finished {m['finished']}, conodb {m['conodb_version']})")
    print(f"Entries: {m['uniprot']['entries']}   previous: {m['previous_release'] or '-'}")
    for stage, info in m["stages"].items():
        print(f"  {stage:<8} {info['finished']}  {info['stats']}")
    print(f"Report: {root / 'releases' / name / 'report.md'}")
    return OK


def _cmd_stages(args) -> int:
    for s in STAGES:
        print(f"{s.name:<8} {s.description}\n         config: {', '.join(s.sections)}\n         outputs: {', '.join(s.outputs)}")
    return OK


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="conodb", description="Build and update the ConoDictor reference database.")
    p.add_argument("--version", action="version", version=f"conodb {__version__}")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("init", help="write an example conodb.toml and labels.tsv")
    s.add_argument("directory", nargs="?", default=".")
    s.set_defaults(func=_cmd_init)

    def with_config(sp):
        sp.add_argument("-c", "--config", type=Path, default=Path("conodb.toml"), help="config file (default: ./conodb.toml)")

    s = sub.add_parser("run", help="(re)build the database; unchanged stages are skipped")
    with_config(s)
    s.add_argument("--force", action="store_true", help="rerun every stage")
    s.add_argument("--from", dest="from_stage", choices=STAGE_NAMES, help="rerun from this stage onwards")
    s.add_argument("--until", choices=STAGE_NAMES, help="stop after this stage")
    s.add_argument("--input", type=Path, help="use a local UniProt JSONL(.gz) instead of downloading")
    s.add_argument("--release", help="release name for --input (default: local-<hash>)")
    s.set_defaults(func=_cmd_run)

    s = sub.add_parser("check", help=f"is a newer UniProt release available? (exit {UPDATE_AVAILABLE} if yes)")
    with_config(s)
    s.set_defaults(func=_cmd_check)

    s = sub.add_parser("status", help="show the latest built release")
    with_config(s)
    s.set_defaults(func=_cmd_status)

    s = sub.add_parser("stages", help="list pipeline stages")
    s.set_defaults(func=_cmd_stages)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    try:
        return args.func(args)
    except CurationNeeded as e:
        logging.error("%s", e)
        return CURATION
    except (FetchError, BuildError, ScoreError, FitError, FileNotFoundError, ValueError) as e:
        logging.error("%s", e)
        return ERROR


if __name__ == "__main__":
    sys.exit(main())
