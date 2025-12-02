"""Command-line interface for sead-sexdiff.

Provides a subcommand to download data from the public SEA-AD S3 bucket
using the `download_data` utility.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

from .data import (
    download_data,
    DEFAULT_BUCKET,
    write_subset_from_folder,
    pseudobulk_and_filter,
)


def _parse_suffixes(values: list[str] | None) -> list[str] | None:
    """Normalize suffix arguments.

    Supports passing multiple `--suffix` flags or a single comma-separated
    value. Returns `None` if no suffixes are provided so that `download_data`
    uses its own default behavior.
    """
    if not values:
        return None
    out: list[str] = []
    for v in values:
        if not v:
            continue
        # allow comma-separated
        parts = [p.strip() for p in v.split(",") if p.strip()]
        out.extend(parts)
    return out or None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sead-sexdiff",
        description="Utilities for SEA-AD sex-difference analyses",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    dl = sub.add_parser(
        "download",
        help="Download data from the SEA-AD public S3 bucket",
        description=(
            "Download public SEA-AD data files via anonymous S3 access. "
            "By default downloads CSV and HTML files under the given prefix. "
            "You may also specify explicit files to download, which overrides suffix filtering."
        ),
    )
    dl.add_argument(
        "--prefix",
        default="PFC/RNAseq/",
        help="S3 prefix (folder path) to search within the bucket (default: %(default)s)",
    )
    dl.add_argument(
        "--bucket",
        default=DEFAULT_BUCKET,
        help="S3 bucket name (default: %(default)s)",
    )
    dl.add_argument(
        "--suffix",
        action="append",
        dest="suffixes",
        metavar="SUFFIX",
        help=(
            "File suffix or pattern to include (e.g., .csv, .html, *.csv). "
            "Repeat the flag or use comma-separated values to provide multiple."
        ),
    )
    dl.add_argument(
        "--file",
        action="append",
        dest="files",
        metavar="FILE",
        help=(
            "Specific file under the prefix to download. Repeat the flag or use "
            "comma-separated values to provide multiple. When provided, this "
            "overrides --suffix filtering. You can pass keys relative to the "
            "prefix (e.g., public_datasets/SEA-AD-cell-annotation.2024-08-27.csv) "
            "or full keys that already start with the prefix."
        ),
    )
    dl.add_argument(
        "--dest-root",
        default="data",
        type=Path,
        help="Local destination root directory (default: %(default)s)",
    )
    dl.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar with percentage and ETA during downloads",
    )

    # subset command: build a concatenated subset and write to disk
    sub_subset = sub.add_parser(
        "subset",
        help="Subset donor .h5ad files by subclass and write a single output",
        description=(
            "Load each donor .h5ad in a folder, subset by obs['Subclass'] == label, "
            "concatenate them, and write one compressed H5AD file."
        ),
    )
    sub_subset.add_argument(
        "--folder",
        default="data/PFC/RNAseq/donor_objects",
        type=Path,
        help=(
            "Input folder that contains donor .h5ad files (default: %(default)s)"
        ),
    )
    sub_subset.add_argument(
        "--subclass-label",
        default="Microglia-PVM",
        help=(
            "Value to match in adata.obs['Subclass'] (default: %(default)s)"
        ),
    )
    sub_subset.add_argument(
        "--out-dir",
        default="results/cleaned_files",
        type=Path,
        help="Directory to write the output file (default: %(default)s)",
    )
    sub_subset.add_argument(
        "--out-name",
        default="microglia_subset.h5ad",
        help="Output filename (default: %(default)s)",
    )
    sub_subset.add_argument(
        "--compression",
        default="gzip",
        help=(
            "Compression for AnnData.write (e.g., gzip). Use 'none' to disable "
            "compression. (default: %(default)s)"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "download":
        suffixes = _parse_suffixes(args.suffixes)
        files = _parse_suffixes(args.files)
        paths: List[Path] = download_data(
            prefix=args.prefix,
            bucket=args.bucket,
            suffixes=suffixes,
            files=files,
            dest_root=args.dest_root,
            show_progress=not args.no_progress,
        )
        for p in paths:
            print(p)
        return 0

    if args.command == "subset":
        compression = None if (str(args.compression).lower() in {"none", "no", "false", "0"}) else args.compression
        out_path = write_subset_from_folder(
            folder=args.folder,
            subclass_label=args.subclass_label,
            out_dir=args.out_dir,
            out_name=args.out_name,
            compression=compression,
        )
        print(out_path)

        # Also compute and save pseudobulk from the produced subset
        try:
            import anndata as ad  # local import to avoid hard dependency in non-subset commands

            subset_adata = ad.read_h5ad(str(out_path))
            pdata = pseudobulk_and_filter(subset_adata)

            pseudobulk_path = Path(args.out_dir) / "microglia_pseudobulk.h5ad"
            pdata.write(str(pseudobulk_path), compression=compression)
            print(pseudobulk_path)
        except Exception as e:
            # Surface a clear error while keeping CLI behavior explicit
            raise RuntimeError(f"Failed to create/save pseudobulk: {e}")
        return 0

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
