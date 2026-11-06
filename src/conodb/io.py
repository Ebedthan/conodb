"""Small I/O helpers: atomic writes, TSV, FASTA, hashing."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any


@contextmanager
def atomic_open(path: Path, mode: str = "w", gz: bool = False):
    """Write to a temp file in the same directory, then rename into place.

    A crashed or interrupted run never leaves a half-written output that a
    later run would mistake for a finished one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    try:
        binary = "b" in mode
        if gz:
            fh = gzip.open(
                tmp, "wb" if binary else "wt", encoding=None if binary else "utf-8"
            )
        elif binary:
            fh = open(tmp, mode)
        else:
            fh = open(tmp, mode, encoding="utf-8", newline="")
        with fh:
            yield fh
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> int:
    n = 0
    with atomic_open(path) as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        w.writeheader()
        for row in rows:
            w.writerow({k: _fmt(row.get(k)) for k in fields})
            n += 1
    return n


def _fmt(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (list, tuple, set)):
        return ";".join(str(x) for x in v)
    return str(v)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> int:
    n = 0
    with atomic_open(path) as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
            n += 1
    return n


def iter_jsonl_gz(path: Path) -> Iterator[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_obj(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, default=str).encode()
    ).hexdigest()


def write_json(path: Path, obj: Any) -> None:
    with atomic_open(path) as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")


def read_json(path: Path) -> Any:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
