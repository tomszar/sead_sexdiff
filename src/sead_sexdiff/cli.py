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
    differential_expression,
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


def _parse_floats_csv(value: str | None) -> list[float] | None:
    """Parse a comma- or space-separated string of numbers into a list of floats.

    Returns None if the input is falsy.
    """
    if not value:
        return None
    parts = []
    # Allow commas and/or whitespace
    for token in value.replace(",", " ").split():
        if not token:
            continue
        parts.append(float(token))
    return parts or None


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
        default=None,
        help=(
            "Value to match in adata.obs['Subclass']. If omitted, a subset file "
            "is produced for EACH Subclass present across the input files."
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
        default=None,
        help=(
            "Output filename (used only when --subclass-label is provided). "
            "If omitted, the filename is inferred from the subclass label "
            "(e.g., 'Microglia-PVM_subset.h5ad'). Ignored when producing all subclasses."
        ),
    )
    sub_subset.add_argument(
        "--compression",
        default="lzf",
        help=(
            "Compression for AnnData.write (e.g., lzf, gzip). Use 'none' to disable "
            "compression. (default: %(default)s)"
        ),
    )

    # differential expression command
    sub_de = sub.add_parser(
        "de",
        help="Run differential expression with pydeseq2 on pseudobulk data",
        description=(
            "Load a pseudobulk .h5ad file and, by default, run DE for ALL supertypes. "
            "Fits a DESeq2 model with design '~Sex + ADNC + Sex:ADNC', computes a "
            "contrast, writes one results table and volcano plot per supertype."
        ),
    )
    sub_de.add_argument(
        "--pseudobulk-file",
        default=Path("results/cleaned_files/microglia_pseudobulk.h5ad"),
        type=Path,
        help="Path to pseudobulk .h5ad file (default: %(default)s)",
    )
    sub_de.add_argument(
        "--supertype",
        default=None,
        help=(
            "Value in obs.Supertype to subset. If omitted, runs DE for ALL "
            "supertypes present in the data."
        ),
    )
    sub_de.add_argument(
        "--results-dir",
        default=Path("results/DE"),
        type=Path,
        help=(
            "Directory to write outputs. A CSV and a volcano PDF are produced per "
            "supertype (default: %(default)s)"
        ),
    )
    sub_de.add_argument(
        "--contrast",
        default=None,
        help=(
            "Contrast vector as comma-separated numbers (e.g., '0,0,0,0,0,0,0,1'). "
            "If omitted, the default contrast (0,0,0,0,0,0,0,1) is used."
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

        print("[subset] Starting subset and concatenate step...")
        print(f"[subset] Folder: {args.folder}")
        print(f"[subset] Subclass label: {args.subclass_label if args.subclass_label is not None else 'ALL (one file per Subclass)'}")
        print(f"[subset] Output directory: {args.out_dir}")
        if args.subclass_label is not None:
            print(f"[subset] Output name: {args.out_name or '(auto from subclass)'}")
        else:
            print("[subset] Output name: (auto: '<Subclass>_subset.h5ad' per subclass)")
        print(f"[subset] Compression: {compression or 'none'}")
        print("[subset] Building concatenated subset (this may take a while)...")
        out_paths = write_subset_from_folder(
            folder=args.folder,
            subclass_label=args.subclass_label,
            out_dir=args.out_dir,
            out_name=args.out_name,
            compression=compression,
        )
        if isinstance(out_paths, list):
            print("[subset] Wrote the following subset files:")
            for p in out_paths:
                print(f"  - {p}")
            print("[subset] Note: Renaming obs column 'Overall AD neuropathological Change' to 'ADNC' is applied before writing when present.")
            print("[subset] Pseudobulk files were generated alongside subset files.")
        else:
            print(f"[subset] Wrote subset file to: {out_paths}")
            print("[subset] Note: Renaming obs column 'Overall AD neuropathological Change' to 'ADNC' is applied before writing when present.")
            print("[subset] Pseudobulk file was generated alongside subset file.")
        return 0

    if args.command == "de":
        print("[de] Running differential expression analysis...")
        contrast_list = _parse_floats_csv(args.contrast)
        if contrast_list is None:
            outputs = differential_expression(
                pseudobulk_file=args.pseudobulk_file,
                supertype=args.supertype,
                results_dir=args.results_dir,
            )
        else:
            outputs = differential_expression(
                pseudobulk_file=args.pseudobulk_file,
                supertype=args.supertype,
                results_dir=args.results_dir,
                contrast=contrast_list,
            )
        if isinstance(outputs, dict):
            # Print a line per supertype
            for st, paths in outputs.items():
                csv_p = paths.get("results_csv")
                pdf_p = paths.get("volcano_pdf")
                if csv_p:
                    print(f"[de] {st} results: {csv_p}")
                if pdf_p:
                    print(f"[de] {st} volcano: {pdf_p}")
        print("[de] Done.")
        return 0

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
