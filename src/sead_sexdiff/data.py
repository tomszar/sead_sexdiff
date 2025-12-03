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
import json

import boto3
from botocore import UNSIGNED
from botocore.client import Config
import pandas as pd
import scanpy as sc
import numpy as np
import random
import decoupler as dc
import anndata as ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


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

    Iterates through files, subsetting by ``Subclass == subclass_label`` and
    additionally excluding samples where ``obs["Overall AD neuropathological Change"] == "Reference"``
    for each file. Concatenates the per-file subsets into a single AnnData while
    aggressively freeing intermediate objects to minimize memory usage.

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
        Concatenation of per-file subsets after applying both filters:
        1) ``obs['Subclass'] == subclass_label`` and
        2) ``obs['Overall AD neuropathological Change'] != 'Reference'``.
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
        # Apply additional subsetting condition
        if "Overall AD neuropathological Change" not in subset.obs:
            try:
                del subset
            finally:
                gc.collect()
            raise KeyError(
                "Column 'Overall AD neuropathological Change' not found in adata.obs"
            )
        subset = subset[
            subset.obs["Overall AD neuropathological Change"] != "Reference"
        ]
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
    subclass_label: str | None = None,
    out_dir: str | Path = "results/cleaned_files",
    out_name: str | None = None,
    compression: str | None = "lzf",
) -> Path | list[Path]:
    """High-level helper: subset donor files and write H5AD outputs.

    Parameters
    ----------
    folder:
        Input folder with donor ``.h5ad`` objects.
    subclass_label:
        Subclass value used for filtering. If ``None`` (default), create one
        subset file per unique value in ``adata.obs['Subclass']`` across the
        input donor files.
    out_dir:
        Output directory where the single concatenated file will be written.
    out_name:
        Output filename. Used only when ``subclass_label`` is provided. If not
        provided, a filename based on the subclass label will be generated
        automatically (e.g., ``Microglia-PVM_subset.h5ad``). Ignored when
        ``subclass_label`` is ``None``.
    compression:
        Compression passed to ``AnnData.write`` (e.g., "gzip").

    Returns
    -------
    pathlib.Path | list[pathlib.Path]
        Path to the written file when a single ``subclass_label`` is provided,
        or a list of paths when processing all subclasses.
    """
    import anndata as ad  # only for typing/IDE hints; actual write uses the object
    from collections import OrderedDict

    def _sanitize(name: str) -> str:
        # Simple filename sanitizer: keep alphanum, dash, underscore; replace others with '_'
        import re

        name = name.strip()
        name = re.sub(r"\s+", "_", name)
        name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        # collapse multiple underscores
        name = re.sub(r"_+", "_", name)
        return name

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # If a specific subclass is requested, process and write a single file
    if subclass_label is not None:
        label = str(subclass_label)
        auto_name = f"{_sanitize(label)}_subset.h5ad"
        out_path = out_dir / (out_name or auto_name)

        adata_concat = subset_and_concat_folder(folder, subclass_label=label)

        # Rename obs column before writing if present
        try:
            if "Overall AD neuropathological Change" in adata_concat.obs.columns:
                adata_concat.obs.rename(
                    columns={"Overall AD neuropathological Change": "ADNC"},
                    inplace=True,
                )
        except Exception:
            # Non-fatal; continue to write
            pass

        adata_concat.write(str(out_path), compression=compression)

        # Immediately create and write pseudobulk, then free memory
        try:
            pdata = pseudobulk_and_filter(adata_concat)
            stem = out_path.stem
            if stem.endswith("_subset"):
                stem = stem[: -len("_subset")]
            pseudobulk_path = out_dir / f"{stem}_pseudobulk.h5ad"
            pdata.write(str(pseudobulk_path), compression=compression)
        finally:
            # Ensure we free large objects aggressively
            try:
                del pdata
            except Exception:
                pass
            try:
                del adata_concat
            except Exception:
                pass
            import gc as _gc
            _gc.collect()
        return out_path

    # Otherwise, discover all subclasses and write one file per subclass
    # Use lightweight reader to avoid loading matrices and cache discovery
    folder = Path(folder)
    files = sorted([p for p in folder.glob("*.h5ad") if p.is_file()])
    if not files:
        raise FileNotFoundError(f"No files matching '*.h5ad' in {folder}")

    # Cache path (under output directory) for subclass discovery
    cache_path = out_dir / ".subclass_index.json"

    def _stat_for(p: Path) -> dict:
        try:
            st = p.stat()
            return {"size": int(st.st_size), "mtime": float(st.st_mtime)}
        except Exception:
            return {"size": None, "mtime": None}

    cached_index: dict[str, dict] | None = None
    if cache_path.exists():
        try:
            cached_index = json.loads(cache_path.read_text())
        except Exception:
            cached_index = None

    # Determine if cache is valid for current files
    cache_valid = False
    if cached_index and isinstance(cached_index, dict):
        files_map = cached_index.get("files") or {}
        try:
            cache_valid = all(
                (
                    str(fp) in files_map
                    and files_map[str(fp)].get("size") == _stat_for(fp)["size"]
                    and abs(files_map[str(fp)].get("mtime", 0.0) - _stat_for(fp)["mtime"]) < 1e-6
                )
                for fp in files
            )
        except Exception:
            cache_valid = False

    # Ordered unique subclasses for deterministic writing order
    subclasses_od: OrderedDict[str, None] = OrderedDict()
    per_file_subclasses: dict[str, list[str]] = {}

    if cache_valid:
        # Rebuild ordered dict by iterating files in current order
        for fp in files:
            vals = cached_index["files"][str(fp)].get("subclasses") or []
            vals = [str(v) for v in vals]
            per_file_subclasses[str(fp)] = vals
            for v in vals:
                subclasses_od.setdefault(v, None)
    else:
        for fp in files:
            try:
                meta = ad.read_h5ad(fp, backed='r')
                obs = meta.obs
                del meta
                vals: list[str] = []
                if obs is not None and "Subclass" in obs.columns:
                    for val in pd.unique(obs["Subclass"].astype(str)):
                        subclasses_od.setdefault(str(val), None)
                        vals.append(str(val))
                per_file_subclasses[str(fp)] = vals
            except Exception:
                # Skip problematic files, continue
                continue

        # Write cache
        try:
            cache_payload = {
                "files": {
                    str(fp): {
                        "size": _stat_for(fp)["size"],
                        "mtime": _stat_for(fp)["mtime"],
                        "subclasses": per_file_subclasses.get(str(fp), []),
                    }
                    for fp in files
                }
            }
            cache_path.write_text(json.dumps(cache_payload, indent=2))
        except Exception:
            # Cache failures are non-fatal
            pass

    if not subclasses_od:
        raise RuntimeError("No subclasses found across input files.")

    written: list[Path] = []
    for label in subclasses_od.keys():
        try:
            print(label)
            adata_concat = subset_and_concat_folder(folder, subclass_label=label)
            try:
                if "Overall AD neuropathological Change" in adata_concat.obs.columns:
                    adata_concat.obs.rename(
                        columns={"Overall AD neuropathological Change": "ADNC"},
                        inplace=True,
                    )
            except Exception:
                pass
            fname = f"{_sanitize(label)}_subset.h5ad"
            out_path = out_dir / fname
            adata_concat.write(str(out_path), compression=compression)
            written.append(out_path)
        except Exception:
            # Continue with next subclass if one fails
            continue

        # Immediately compute and write pseudobulk, then free memory
        try:
            pdata = pseudobulk_and_filter(adata_concat)
            stem = out_path.stem
            if stem.endswith("_subset"):
                stem = stem[: -len("_subset")]
            pseudobulk_path = out_dir / f"{stem}_pseudobulk.h5ad"
            pdata.write(str(pseudobulk_path), compression=compression)
        finally:
            # Ensure we free large objects aggressively
            try:
                del pdata
            except Exception:
                pass
            try:
                del adata_concat
            except Exception:
                pass
            import gc as _gc
            _gc.collect()

    if not written:
        raise RuntimeError("Failed to create any subclass subset files.")
    return written



def prep_anndata(adata_, *, min_cells: int = 3, min_genes: int = 200):
    """Prepare an AnnData object exported from R/Seurat for Scanpy workflows.

    This helper fixes common conversion issues when moving between R and Python
    and applies basic QC filtering:
    - Ensure ``obs`` column dtypes are properly interpreted (e.g. byte strings
      become UTF-8 strings) by reconstructing the object through a temporary
      ``DataFrame`` merge.
    - Filter out genes that are expressed in a small number of cells and cells
      with too few detected genes.

    Parameters
    ----------
    adata_ : anndata.AnnData
        Input AnnData object. Typically created by reading a ``.h5ad`` that was
        exported from R/Seurat or other pipelines where categorical/text columns
        may arrive as bytes or mixed dtypes.
    min_cells : int, default 3
        Keep genes that are expressed in at least ``min_cells`` cells. Passed to
        ``scanpy.pp.filter_genes`` as ``min_cells``.
    min_genes : int, default 200
        Keep cells that have at least ``min_genes`` detected genes. Passed to
        ``scanpy.pp.filter_cells`` as ``min_genes``.

    Returns
    -------
    anndata.AnnData
        A new AnnData object with fixed dtypes and filtered cells/genes.

    Notes
    -----
    - The function returns a fresh AnnData instance; the input object is not
      modified in place.
    - The dtype-fix step reconstructs ``X``, ``obs`` and ``var`` via a
      ``pandas.DataFrame`` to coerce potential byte-string columns to proper
      Python strings and to ensure consistent dtypes across columns.
    """

    def fix_dtypes(adata_):
        df = pd.DataFrame(adata_.X.A, index=adata_.obs_names, columns=adata_.var_names)
        df = df.join(adata_.obs)
        return sc.AnnData(df[adata_.var_names], obs=df.drop(columns=adata_.var_names))

    adata_ = fix_dtypes(adata_)
    sc.pp.filter_genes(adata_, min_cells=min_cells)
    sc.pp.filter_cells(adata_, min_genes=min_genes)
    return adata_


def pseudobulk_and_filter(
    adata,
    *,
    sample_col: str = "Donor ID",
    groups_col: str = "Supertype",
    mode: str = "sum",
    layer: str | None = "UMIs",
    min_cells: int = 10,
    min_counts: int = 1000,
    pseudobulk_kwargs: dict | None = None,
    filter_kwargs: dict | None = None,
):
    """Create pseudobulk profiles and filter low-quality samples using decoupler.

    This helper wraps two decoupler preprocessing utilities in sequence:
    1) ``dc.pp.pseudobulk`` to aggregate single-cell counts into sample/group
       pseudobulk profiles.
    2) ``dc.pp.filter_samples`` to filter out samples with too few cells or
       total counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input single-cell AnnData object.
    sample_col : str, default "Donor ID"
        Column in ``adata.obs`` identifying the sample (e.g., donor ID).
    groups_col : str, default "Supertype"
        Column in ``adata.obs`` defining the biological group/cluster to
        aggregate within (e.g., cell type).
    mode : str, default "sum"
        Aggregation mode for ``dc.pp.pseudobulk`` (e.g., "sum", "mean").
    layer : str | None, default "UMIs"
        Name of the layer to use for counts. If ``None``, use ``adata.X``.
    min_cells : int, default 10
        Minimum number of cells required per sample/group to keep in
        ``dc.pp.filter_samples``.
    min_counts : int, default 1000
        Minimum total counts required per sample/group to keep in
        ``dc.pp.filter_samples``.
    pseudobulk_kwargs : dict | None, optional
        Extra keyword arguments forwarded to ``dc.pp.pseudobulk``.
    filter_kwargs : dict | None, optional
        Extra keyword arguments forwarded to ``dc.pp.filter_samples``.

    Returns
    -------
    anndata.AnnData
        The pseudobulk AnnData after filtering.

    Notes
    -----
    - Requires the ``decoupler`` package to be installed and importable as
      ``decoupler``.
    """

    pb_kwargs = dict(
        adata=adata,
        sample_col=sample_col,
        groups_col=groups_col,
        mode=mode,
        layer=layer,
    )
    if pseudobulk_kwargs:
        pb_kwargs.update(pseudobulk_kwargs)

    pdata = dc.pp.pseudobulk(**pb_kwargs)

    fs_kwargs = dict(min_cells=min_cells, min_counts=min_counts)
    if filter_kwargs:
        fs_kwargs.update(filter_kwargs)

    dc.pp.filter_samples(pdata, **fs_kwargs)
    return pdata


def differential_expression(
    pseudobulk_file: str | Path = "results/cleaned_files/microglia_pseudobulk.h5ad",
    *,
    supertype: str | None = None,
    results_dir: str | Path = "results/DE",
    contrast: np.ndarray | list[float] | tuple[float, ...] = (
        0, 0, 0, 0, 0, 0, 0, 1
    ),
):
    """Run DESeq2 differential expression on pseudobulked data.

    Steps:
    - Load pseudobulk AnnData file.
    - For each category in ``obs.Supertype`` (or a specific one if provided),
      create a ``DeseqDataSet`` with design ``~Sex + ADNC + Sex:ADNC``.
    - Fit the model, compute statistics for a single contrast, and summarize.
    - Save a results table and a volcano plot per supertype.

    Parameters
    ----------
    pseudobulk_file : str | Path
        Path to the pseudobulk AnnData file (``.h5ad``).
    supertype : str | None
        If provided, limit the analysis to this value in ``obs.Supertype``.
        If ``None``, run for all supertypes present in the data.
    results_dir : str | Path, default "results/DE"
        Directory to write results into. A CSV and a PDF are produced per
        supertype as ``<Supertype>_deseq2_results.csv`` and
        ``<Supertype>_volcano.pdf``.
    contrast : array-like of float, default (0,0,0,0,0,0,0,1)
        Contrast vector passed to ``pydeseq2.DeseqStats``. Provide a sequence
        of numbers matching your design matrix. Defaults to the requested
        vector with the interaction term.

    Behavior on non full-rank design
    --------------------------------
    If ``pydeseq2`` emits a ``UserWarning`` stating that "The design matrix is
    not full rank, so the model cannot be fitted", this function skips that
    supertype and proceeds to the next one without writing outputs for it.

    Returns
    -------
    dict[str, dict[str, Path]]
        Mapping of ``supertype`` to a dict with keys ``results_csv`` and
        ``volcano_pdf``.
    """

    pseudobulk_file = Path(pseudobulk_file)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load pseudobulk AnnData
    pdata = ad.read_h5ad(str(pseudobulk_file))

    # Validate columns
    if "Supertype" not in pdata.obs.columns:
        raise KeyError("Column 'Supertype' not found in pseudobulk.obs")
    for col in ("Sex", "ADNC"):
        if col not in pdata.obs.columns:
            raise KeyError(f"Column '{col}' not found in pseudobulk.obs")

    # Determine supertypes to process
    if supertype is None:
        supertypes = pd.Index(sorted(pd.unique(pdata.obs["Supertype"].astype(str))))
    else:
        supertypes = pd.Index([supertype])

    outputs: dict[str, dict[str, Path]] = {}

    for st in supertypes:
        st_mask = pdata.obs["Supertype"].astype(str) == st
        micro_pdata = pdata[st_mask].copy()
        if micro_pdata.n_obs == 0:
            # Skip empty groups (in case of mismatched dtype)
            continue

        # Ensure categorical covariates are strings/categorical
        for col in ("Sex", "ADNC"):
            # Convert bytes to str if any
            if micro_pdata.obs[col].dtype == object:
                micro_pdata.obs[col] = micro_pdata.obs[col].apply(
                    lambda x: x.decode() if isinstance(x, (bytes, bytearray)) else x
                )
            micro_pdata.obs[col] = micro_pdata.obs[col].astype("category")

        # Create DESeq2 dataset and fit
        # If pydeseq2 emits a UserWarning that the design matrix is not full rank,
        # skip this supertype and continue with the next one.
        import warnings

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "error",
                    message=r"The design matrix is not full rank, so the model cannot be fitted.*",
                    category=UserWarning,
                )
                dds = DeseqDataSet(adata=micro_pdata, design="~Sex + ADNC + Sex:ADNC")
                # Some versions could defer the check to deseq2(); keep the same guard
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "error",
                        message=r"The design matrix is not full rank, so the model cannot be fitted.*",
                        category=UserWarning,
                    )
                    dds.deseq2()
        except UserWarning as uw:  # design matrix not full rank
            print(
                f"[de] Skipping supertype '{st}' due to non-full-rank design: {uw}")
            continue

        # Use provided contrast vector
        contrast_vec = np.array(contrast, dtype=float)
        ds = DeseqStats(dds, contrast=contrast_vec)
        ds.summary()

        # Save results table
        safe_super = str(st).replace("/", "-")
        out_csv = results_dir / f"{safe_super}_deseq2_results.csv"
        ds.results_df.to_csv(out_csv, index=True)

        # Volcano plot per supertype with name including supertype
        volcano_pdf = results_dir / f"{safe_super}_volcano.pdf"
        dc.pl.volcano(ds.results_df, x="log2FoldChange", y="pvalue", save=str(volcano_pdf))

        outputs[str(st)] = {"results_csv": out_csv, "volcano_pdf": volcano_pdf}

    return outputs

__all__ = [
    "DEFAULT_BUCKET",
    "download_data",
    "subset_adata",
    "subset_and_concat_folder",
    "write_subset_from_folder",
    "prep_anndata",
    "pseudobulk_and_filter",
    "differential_expression",
]