"""UniProt REST client: release probe and paginated download.

Uses only the standard library. Pagination follows the `Link: <...>;
rel="next"` header (UniProt's cursor pagination), retries transient errors
with exponential backoff, and checks the number of entries received
against `X-Total-Results` so a truncated download can never be used.
"""

from __future__ import annotations

import json
import logging
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from . import __version__
from .config import SourceConfig
from .io import atomic_open

log = logging.getLogger(__name__)

USER_AGENT = f"conodb/{__version__} (ConoDictor database builder; python-urllib)"
_RETRY_STATUS = {429, 500, 502, 503, 504}
_LINK_NEXT = re.compile(r'<([^>]+)>\s*;\s*rel="next"')


class FetchError(RuntimeError):
    pass


@dataclass
class ReleaseInfo:
    release: str
    release_date: str
    total: int

    @property
    def dirname(self) -> str:
        return self.release


def _url(cfg: SourceConfig, size: int) -> str:
    params = {"query": cfg.query, "format": "json", "size": str(size)}
    return f"{cfg.base_url}?{urllib.parse.urlencode(params)}"


def _get(url: str, cfg: SourceConfig) -> tuple[bytes, dict[str, str]]:
    delay = 2.0
    for attempt in range(1, cfg.retries + 1):
        req = urllib.request.Request(
            url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"}
        )
        try:
            with urllib.request.urlopen(req, timeout=cfg.timeout) as resp:
                return resp.read(), {k.lower(): v for k, v in resp.headers.items()}
        except urllib.error.HTTPError as e:
            if e.code not in _RETRY_STATUS or attempt == cfg.retries:
                body = e.read()[:500].decode("utf-8", "replace")
                raise FetchError(f"HTTP {e.code} for {url}: {body}") from e
            wait = float(e.headers.get("Retry-After", delay)) if e.headers else delay
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt == cfg.retries:
                raise FetchError(f"Network error for {url}: {e}") from e
            wait = delay
        log.warning(
            "Request failed (attempt %d/%d), retrying in %.0fs",
            attempt,
            cfg.retries,
            wait,
        )
        time.sleep(wait)
        delay = min(delay * 2, 120)
    raise FetchError(f"Giving up on {url}")  # pragma: no cover


def next_link(headers: dict[str, str]) -> str | None:
    m = _LINK_NEXT.search(headers.get("link", ""))
    return m.group(1) if m else None


def _release(headers: dict[str, str], total: int) -> ReleaseInfo:
    release = headers.get("x-uniprot-release") or f"snapshot-{date.today().isoformat()}"
    return ReleaseInfo(
        release=release.replace("/", "_"),
        release_date=headers.get("x-uniprot-release-date", ""),
        total=total,
    )


def probe(cfg: SourceConfig) -> ReleaseInfo:
    """One tiny request: current UniProt release and number of matching entries."""
    _, headers = _get(_url(cfg, 1), cfg)
    total = int(headers.get("x-total-results", "-1"))
    info = _release(headers, total)
    log.info(
        "UniProt release %s (%s): %d entries match the query",
        info.release,
        info.release_date or "?",
        total,
    )
    return info


def fetch_all(cfg: SourceConfig, out: Path) -> ReleaseInfo:
    """Download every matching entry as JSON lines (gzip) and return release info."""
    url: str | None = _url(cfg, cfg.page_size)
    n = 0
    info: ReleaseInfo | None = None
    with atomic_open(out, gz=True) as fh:
        while url:
            body, headers = _get(url, cfg)
            if info is None:
                info = _release(headers, int(headers.get("x-total-results", "-1")))
            elif (
                headers.get("x-uniprot-release", info.release).replace("/", "_")
                != info.release
            ):
                raise FetchError("UniProt release changed during download; rerun later")
            for entry in json.loads(body).get("results", []):
                fh.write(json.dumps(entry, separators=(",", ":")) + "\n")
                n += 1
            log.info(
                "  downloaded %d / %s entries",
                n,
                info.total if info.total >= 0 else "?",
            )
            url = next_link(headers)
        assert info is not None
        if info.total >= 0 and n != info.total:
            raise FetchError(
                f"Expected {info.total} entries, received {n}; download incomplete"
            )
    info.total = n
    return info
