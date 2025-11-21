"""Data utilities for downloading resources from the SEA-AD S3 bucket.

This module provides a convenience function `download_data()` to fetch public
files from the SEA-AD single-cell profiling S3 bucket without requiring AWS
credentials. It supports filtering by suffix (e.g., ``.csv``, ``.html``), and
will always attempt to fetch the ``Changelog.html`` file within the specified
prefix.

Bucket browser: https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html

Example
-------
>>> from sead_sexdiff.data import download_data
>>> downloaded = download_data(prefix="PFC/RNAseq/", suffixes=[".csv", ".html"])  # doctest: +SKIP
>>> print(downloaded)  # doctest: +SKIP

Notes
-----
- The bucket is public and does not require credentials. We configure boto3 to
  use unsigned requests via botocore.
- Files are downloaded preserving their S3 key structure under the local
  destination directory (default: ``data/``).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Iterable, List

import boto3
from botocore import UNSIGNED
from botocore.client import Config


DEFAULT_BUCKET = "sea-ad-single-cell-profiling"


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _matches_suffixes(key: str, suffixes: Iterable[str] | None) -> bool:
    if not suffixes:
        return True
    for suf in suffixes:
        # Normalize: allow users to pass extensions with or without leading dot
        if not suf:
            # empty string matches everything (but ignore for safety)
            continue
        if suf.startswith("*"):
            # Simple wildcard pattern like *.csv — fallback to endswith after '*'
            suf = suf.lstrip("*")
        if not suf.startswith(".") and "/" not in suf:
            # If it's likely an extension without dot, add dot for convenience
            candidate = "." + suf
        else:
            candidate = suf
        if key.endswith(candidate) or key == suf:
            return True
    return False


def download_data(
    prefix: str = "PFC/RNAseq/",
    *,
    bucket: str = DEFAULT_BUCKET,
    suffixes: Iterable[str] | None = (".csv", ".html"),
    files: Iterable[str] | None = None,
    dest_root: str | Path = "data",
    show_progress: bool = True,
) -> List[Path]:
    """Download public SEA-AD data from S3 via boto3.

    Parameters
    ----------
    prefix:
        S3 prefix (folder-like path) to search within the bucket. For example
        "PFC/RNAseq/".
    bucket:
        S3 bucket name. Defaults to "sea-ad-single-cell-profiling".
    suffixes:
        File suffixes/extensions to download (e.g., [".csv", ".html"]). If
        ``None`` or empty, all objects under the prefix will be downloaded.
        You can also pass patterns like "*.csv" or bare extensions like
        "csv". Matching is performed by simple endswith checks.
    files:
        Specific file paths under the given ``prefix`` to download. These may
        be provided either as keys relative to ``prefix`` (e.g.,
        "public_datasets/SEA-AD-cell-annotation.2024-08-27.csv") or as full
        S3 keys starting with the same ``prefix``. When this option is
        provided and non-empty, it overrides ``suffixes`` filtering and only
        the specified files (plus the ``Changelog.html``, if present) are
        attempted.
    dest_root:
        Local root directory where files will be saved. The S3 key structure is
        preserved under this directory.

    Returns
    -------
    List[pathlib.Path]
        A list of local file paths that were downloaded.

    Notes
    -----
    This function uses unsigned S3 requests, which is appropriate for public
    buckets. If the bucket ever becomes private, authentication would be
    required and this approach would need updating.
    """

    # Normalize inputs
    if prefix and not prefix.endswith("/"):
        # Keep behavior consistent with folder-like prefixes
        prefix = prefix + "/"

    dest_root = Path(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)

    # Create an anonymous S3 client (no credentials needed for public buckets)
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    downloaded: List[Path] = []

    # Internal progress tracker
    class _Progress:
        def __init__(self, total: int | None):
            self.total = total if (isinstance(total, int) and total > 0) else None
            self.seen = 0
            self.start = time.time()
            self._last_print = 0.0

        def cb(self, bytes_amount: int):
            # Callback invoked by boto3 on each chunk
            self.seen += int(bytes_amount)
            if not show_progress:
                return
            now = time.time()
            # Throttle updates to ~10 Hz
            if now - self._last_print < 0.1:
                return
            self._last_print = now
            self._print_line(end="\r")

        def _eta_str(self, rate: float) -> str:
            if rate <= 0:
                return "--:--"
            remaining = (self.total - self.seen) / rate if self.total else 0
            m, s = divmod(int(max(0, remaining)), 60)
            if m >= 60:
                h, m = divmod(m, 60)
                return f"{h:d}h{m:02d}m"
            return f"{m:02d}:{s:02d}"

        def _human(self, n: int) -> str:
            # Simple human-readable bytes
            units = ["B", "KB", "MB", "GB", "TB"]
            size = float(n)
            for u in units:
                if size < 1024 or u == units[-1]:
                    return f"{size:.1f}{u}"
                size /= 1024
            return f"{n}B"

        def _print_line(self, end: str = "\n"):
            elapsed = max(1e-6, time.time() - self.start)
            rate = self.seen / elapsed
            rate_str = self._human(int(rate)) + "/s"
            if self.total:
                pct = (self.seen / self.total) * 100
                eta = self._eta_str(rate)
                line = (
                    f"Downloading: {self._human(self.seen)} / {self._human(self.total)} "
                    f"({pct:5.1f}%) | ETA {eta} | {rate_str}"
                )
            else:
                line = f"Downloading: {self._human(self.seen)} | {rate_str}"
            sys.stdout.write(line + end)
            sys.stdout.flush()

        def finish(self):
            if show_progress:
                # Ensure a final full line with 100% if total known
                if self.total:
                    self.seen = max(self.seen, self.total)
                self._print_line(end="\n")

    # Build a queue of downloads with known sizes when possible
    download_queue: list[tuple[str, Path, int | None]] = []  # (key, local_path, size)

    # Always consider Changelog.html
    changelog_key = f"{prefix}Changelog.html"
    try:
        local_changelog = dest_root / changelog_key
        remote_size: int | None = None
        try:
            head = s3.head_object(Bucket=bucket, Key=changelog_key)
            remote_size = int(head.get("ContentLength", -1))
            if remote_size < 0:
                remote_size = None
        except Exception:
            remote_size = None

        # Skip if up-to-date
        try:
            if local_changelog.exists() and remote_size is not None:
                if local_changelog.stat().st_size == remote_size:
                    pass  # do not enqueue
                else:
                    download_queue.append((changelog_key, local_changelog, remote_size))
            else:
                download_queue.append((changelog_key, local_changelog, remote_size))
        except Exception:
            download_queue.append((changelog_key, local_changelog, remote_size))
    except Exception:
        # ignore issues with changelog discovery
        pass

    # If specific files are requested, prepare those and process early return
    if files:
        for f in files:
            if not f:
                continue
            rel = str(f).lstrip("/")
            key = rel if rel.startswith(prefix) else f"{prefix}{rel}"
            local_path = dest_root / key

            remote_size: int | None = None
            try:
                head = s3.head_object(Bucket=bucket, Key=key)
                rs = int(head.get("ContentLength", -1))
                remote_size = rs if rs >= 0 else None
            except Exception:
                remote_size = None

            # Skip re-download if size matches
            try:
                if local_path.exists() and remote_size is not None:
                    if local_path.stat().st_size == remote_size:
                        continue
            except Exception:
                pass

            download_queue.append((key, local_path, remote_size))

        # Execute downloads with progress
        total_bytes = sum(sz for _, _, sz in download_queue if isinstance(sz, int)) or None
        prog = _Progress(total_bytes)
        for key, local_path, _sz in download_queue:
            try:
                _ensure_parent_dir(local_path)
                s3.download_file(bucket, key, str(local_path), Callback=prog.cb)
                downloaded.append(local_path)
            except Exception:
                continue
        prog.finish()
        return downloaded

    # List all objects under the prefix and download those matching the suffixes
    paginator = s3.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=bucket, Prefix=prefix)

    for page in pages:
        contents = page.get("Contents") or []
        for obj in contents:
            key = obj.get("Key")
            if not key or key.endswith("/"):
                continue
            if not _matches_suffixes(key, suffixes):
                continue
            local_path = dest_root / key
            # Decide whether to download and capture size
            remote_size = obj.get("Size")
            try:
                if local_path.exists() and isinstance(remote_size, int):
                    if local_path.stat().st_size == remote_size:
                        continue
            except OSError:
                pass
            download_queue.append((key, local_path, remote_size if isinstance(remote_size, int) else None))

    # Execute queued downloads with progress
    total_bytes = sum(sz for _, _, sz in download_queue if isinstance(sz, int)) or None
    prog = _Progress(total_bytes)
    for key, local_path, _sz in download_queue:
        try:
            _ensure_parent_dir(local_path)
            s3.download_file(bucket, key, str(local_path), Callback=prog.cb)
            downloaded.append(local_path)
        except Exception:
            continue
    prog.finish()

    return downloaded


__all__ = ["download_data"]
