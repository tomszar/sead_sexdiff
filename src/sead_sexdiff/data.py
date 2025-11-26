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
    - To download all files under a folder/prefix (e.g. the donor_objects
      folder), pass ``suffixes=None``. For example:
      ``download_data(prefix="PFC/RNAseq/donor_objects/", suffixes=None)``.
    - This function uses unsigned S3 requests, which is appropriate for public
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


def subset_adata(file_path: str | Path, *, subclass_label: str = "Microglia-PVM"):
    """Read an AnnData H5AD file and return a subset by ``obs['Subclass']``.

    Parameters
    ----------
    file_path:
        Path to the ``.h5ad`` file to read.
    subclass_label:
        The label to match in ``adata.obs['Subclass']``. Defaults to
        "Microglia-PVM".

    Returns
    -------
    anndata.AnnData
        The subset ``AnnData`` object. The original loaded object is deleted
        to free memory before returning.
    """
    # Import locally to avoid heavy import on module import
    import gc
    import anndata as ad

    file_path = Path(file_path)
    adata = ad.read_h5ad(str(file_path))
    if "Subclass" not in adata.obs:
        try:
            del adata
        finally:
            gc.collect()
        raise KeyError("Column 'Subclass' not found in adata.obs")

    subset = adata[adata.obs["Subclass"] == subclass_label].copy()
    # Drop the original to free memory
    try:
        del adata
    finally:
        gc.collect()
    return subset


def subset_and_concat_folder(
    folder: str | Path,
    *,
    subclass_label: str = "Microglia-PVM",
    pattern: str = "*.h5ad",
    recursive: bool = False,
):
    """Create a concatenated subset from all H5AD files in a folder.

    Iterates through files, subsetting by ``Subclass == subclass_label`` for
    each file, concatenating them into a single AnnData while aggressively
    freeing intermediate objects to minimize memory usage.

    Parameters
    ----------
    folder:
        Directory containing input ``.h5ad`` files (e.g.
        ``data/PFC/RNAseq/donor_objects``).
    subclass_label:
        Label to match in ``obs['Subclass']`` during subsetting.
    pattern:
        Glob pattern of files to include (default: ``*.h5ad``).
    recursive:
        If True, search subdirectories recursively using ``rglob``.

    Returns
    -------
    anndata.AnnData
        Concatenation of per-file subsets.
    """
    import gc
    from typing import Iterable as _Iter
    import anndata as ad

    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder}")

    files: _Iter[Path]
    files = folder.rglob(pattern) if recursive else folder.glob(pattern)
    files = sorted([p for p in files if p.is_file()])
    if not files:
        raise FileNotFoundError(f"No files matching {pattern!r} in {folder}")

    result: ad.AnnData | None = None
    for fp in files:
        subset = subset_adata(fp, subclass_label=subclass_label)
        if result is None:
            result = subset
        else:
            # Concatenate incrementally, then drop intermediates
            new_result = ad.concat([result, subset], join="outer", merge="unique")
            try:
                del result
                del subset
            finally:
                gc.collect()
            result = new_result

    assert result is not None  # for type checkers
    return result


def write_subset_from_folder(
    folder: str | Path = "data/PFC/RNAseq/donor_objects",
    *,
    subclass_label: str = "Microglia-PVM",
    out_dir: str | Path = "results/cleaned_files",
    out_name: str = "microglia_subset.h5ad",
    compression: str | None = "gzip",
) -> Path:
    """High-level helper: subset all files in ``folder`` and write one H5AD.

    Parameters
    ----------
    folder:
        Input folder with donor ``.h5ad`` objects.
    subclass_label:
        Subclass value used for filtering; defaults to "Microglia-PVM".
    out_dir:
        Output directory where the single concatenated file will be written.
    out_name:
        Filename for the output (default: "microglia_subset.h5ad").
    compression:
        Compression passed to ``AnnData.write`` (e.g., "gzip").

    Returns
    -------
    pathlib.Path
        The path to the written file.
    """
    import anndata as ad  # only for typing/IDE hints; actual write uses the object

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / out_name

    adata_concat = subset_and_concat_folder(folder, subclass_label=subclass_label)
    adata_concat.write(str(out_path), compression=compression)
    return out_path


__all__ = [
    "DEFAULT_BUCKET",
    "download_data",
    "subset_adata",
    "subset_and_concat_folder",
    "write_subset_from_folder",
]
