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

    # Always attempt to download Changelog.html in the given prefix
    changelog_key = f"{prefix}Changelog.html"
    try:
        local_changelog = dest_root / changelog_key
        _ensure_parent_dir(local_changelog)
        s3.download_file(bucket, changelog_key, str(local_changelog))
        downloaded.append(local_changelog)
    except Exception:
        # Non-fatal if not present or any transient error
        pass

    # If specific files are requested, download those and return early
    if files:
        for f in files:
            if not f:
                continue
            # Normalize input to an S3 key under the prefix
            rel = str(f).lstrip("/")
            key = rel if rel.startswith(prefix) else f"{prefix}{rel}"

            local_path = dest_root / key
            # Attempt to skip re-download when sizes match by peeking head
            try:
                if local_path.exists():
                    try:
                        head = s3.head_object(Bucket=bucket, Key=key)
                        remote_size = int(head.get("ContentLength", -1))
                        if remote_size > -1 and local_path.stat().st_size == remote_size:
                            # already up-to-date
                            continue
                    except Exception:
                        # If HEAD fails for any reason, fall back to download
                        pass

                _ensure_parent_dir(local_path)
                s3.download_file(bucket, key, str(local_path))
                downloaded.append(local_path)
            except Exception:
                # Non-fatal; skip missing or inaccessible keys
                continue

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
            # Skip re-download if file exists and size matches
            try:
                if local_path.exists():
                    # Compare sizes when available
                    local_size = local_path.stat().st_size
                    remote_size = obj.get("Size")
                    if isinstance(remote_size, int) and remote_size == local_size:
                        continue
            except OSError:
                # If any filesystem error occurs, attempt to re-download
                pass

            _ensure_parent_dir(local_path)
            s3.download_file(bucket, key, str(local_path))
            downloaded.append(local_path)

    return downloaded


__all__ = ["download_data"]
