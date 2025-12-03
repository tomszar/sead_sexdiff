# sead_sexdiff

Utilities and CLI for SEA-AD sex-difference analyses using Python >= 3.12. Includes tools to download data from the public SEA-AD S3 bucket and to subset/concatenate AnnData (.h5ad) donor files.

Project layout:

```
.
├── data/            # place raw/input data here (ignored by package build)
├── results/         # analysis outputs, figures, tables, etc.
├── src/
│   └── sead_sexdiff/
│       └── __init__.py
├── pyproject.toml   # project configuration (uv/pip compatible)
├── LICENSE
└── README.md
```

Quick start (with uv):

1) Ensure you have Python 3.12+ installed.

2) Install uv (if not already):

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

3) Create and use a virtual environment and install deps:

```
uv venv --python 3.12
source .venv/bin/activate  # or .venv\Scripts\activate on Windows
uv pip install -e .
```

This installs the project dependencies declared in `pyproject.toml` (e.g., `pandas`, `numpy`, `boto3/botocore`, `anndata`, `scanpy`).

### How to use

After installing the package (see Quick start), you can download data with the CLI:

- Default behavior (downloads CSV and HTML under the default prefix `PFC/RNAseq/` into `data/`):

```
sead-sexdiff download
```

- Specify a prefix and multiple suffixes (repeatable and/or comma-separated):

```
sead-sexdiff download --prefix PFC/RNAseq/ --suffix .csv --suffix .html
# or
sead-sexdiff download --suffix .csv,.html
```

- Choose explicit files to download (overrides suffix filtering). You can pass
  keys relative to the prefix or full keys that already start with the prefix.

```
sead-sexdiff download \
  --prefix PFC/RNAseq/ \
  --file public_datasets/SEA-AD-cell-annotation.2024-08-27.csv \
  --file SEAAD_A9_RNAseq_all-nuclei_metadata.2024-02-13.csv

# or comma-separated
sead-sexdiff download --file public_datasets/SEA-AD-cell-annotation.2024-08-27.csv,SEAAD_A9_RNAseq_all-nuclei_metadata.2024-02-13.csv
```

- Change destination directory and bucket if needed:

```
sead-sexdiff download --dest-root ./data --bucket sea-ad-single-cell-profiling
```

The command prints the local paths of files it downloaded. It uses anonymous S3 access for the public SEA-AD bucket and always attempts to fetch `Changelog.html` under the chosen prefix.

### New: Subset donor AnnData files and write pseudobulk(s) (CLI)

Use the `subset` command to load each donor `.h5ad` in a folder, filter by subclass, apply rudimentary QC, concatenate, and write outputs to disk.

- Single-subclass mode (when `--subclass-label` is provided):
  - Concatenated subset: `results/cleaned_files/<Subclass>_subset.h5ad` by default, or the value of `--out-name` if given.
  - Pseudobulk AnnData is also written automatically as `results/cleaned_files/<Subclass>_pseudobulk.h5ad`.

- ALL-subclasses mode (default when `--subclass-label` is omitted):
  - One concatenated subset file per unique `obs['Subclass']`, named `results/cleaned_files/<Subclass>_subset.h5ad`.
  - One pseudobulk per subset is also generated, named `results/cleaned_files/<Subclass>_pseudobulk.h5ad`.

```
sead-sexdiff subset \
  --folder data/PFC/RNAseq/donor_objects \
  --subclass-label Microglia-PVM \
  --out-dir results/cleaned_files \
  --out-name microglia_subset.h5ad \
  --compression gzip
```

- Disable compression if desired: `sead-sexdiff subset --compression none`
- Memory: files are processed one-by-one; intermediates are freed to keep memory usage low.
- QC: when reading each donor object, two basic filters are applied before concatenation:
  - filter cells with too few detected genes: `sc.pp.filter_cells(adata, min_genes=200)`
  - filter genes detected in very few cells: `sc.pp.filter_genes(adata, min_cells=3)`
- Pseudobulk details: created via `decoupler` using defaults (`sample_col='Donor ID'`, `groups_col='Supertype'`, `mode='sum'`, `layer='UMIs'`, `min_cells=10`, `min_counts=1000`). If this step fails (e.g., missing dependency), the CLI raises a clear error. In ALL-subclasses mode, pseudobulk computation runs once per produced subset file.

### New: Differential expression (CLI)

Run DESeq2 differential expression on the pseudobulk file using design `~Sex + ADNC + Sex:ADNC`.

- By default, the command runs DE for ALL supertypes found in `obs.Supertype`.
- Outputs are written to `results/DE`, one CSV and one volcano plot per supertype.

Basic usage (runs all supertypes):

```
sead-sexdiff de
```

Options:

- `--pseudobulk-file`: path to the pseudobulk `.h5ad` (default: `results/cleaned_files/microglia_pseudobulk.h5ad`)
- `--supertype`: value in `obs.Supertype` to subset. If omitted, runs DE for ALL supertypes.
- `--results-dir`: output directory for per-supertype results (default: `results/DE`)
- `--contrast`: contrast vector as comma- or space-separated numbers. If omitted, the default `(0,0,0,0,0,0,0,1)` is used.

Examples:

- Run all supertypes to `results/DE`:

```
sead-sexdiff de
```

- Limit to a single Supertype and customize contrast:

```
sead-sexdiff de \
  --pseudobulk-file results/cleaned_files/microglia_pseudobulk.h5ad \
  --supertype Micro-PVM_2 \
  --results-dir results/DE \
  --contrast 0,0,0,0,0,0,0,1
```

Outputs (per supertype, with `/` in names replaced by `-`):

- Results table CSV: `results/DE/<Supertype>_deseq2_results.csv`
- Volcano plot PDF: `results/DE/<Supertype>_volcano.pdf`

Notes:

- If the DESeq2 design matrix is not full rank for a given supertype (pydeseq2
  emits a UserWarning that the model cannot be fitted), that supertype is
  skipped and no outputs are written for it. The analysis proceeds with the
  remaining supertypes.

### Python API examples

In addition to the CLI, you can call the underlying functions from Python.

- Download all files under a specific folder/prefix (e.g., `donor_objects`). Pass `suffixes=None` to download everything:

```python
from sead_sexdiff.data import download_data

downloaded_paths = download_data(
    prefix="PFC/RNAseq/donor_objects/",
    suffixes=None,  # None => download all files under the prefix
)
```

- Subset a single donor `.h5ad` by subclass and free the original object:

```python
from sead_sexdiff.data import subset_adata

adata_sub = subset_adata(
    "data/PFC/RNAseq/donor_objects/H19.33.004_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad",
    subclass_label="Microglia-PVM",
)
```

Note: `subset_adata` applies the same rudimentary QC by default (cells `min_genes=200`, genes `min_cells=3`).

- Subset every donor file in a folder and concatenate them:

```python
from sead_sexdiff.data import subset_and_concat_folder

adata_concat = subset_and_concat_folder(
    "data/PFC/RNAseq/donor_objects",
    subclass_label="Microglia-PVM",
)
```

- Build and write the concatenated subset to disk (gzip by default):

```python
from sead_sexdiff.data import write_subset_from_folder

out_path = write_subset_from_folder(
    folder="data/PFC/RNAseq/donor_objects",
    subclass_label="Microglia-PVM",
    out_dir="results/cleaned_files",
    out_name="microglia_subset.h5ad",
    compression="gzip",
)
print(out_path)
```

Notes
- The subsetting functions expect an `obs['Subclass']` column in each AnnData object. An error is raised if it is missing.
- To download “all files” under a prefix from Python, use `download_data(..., suffixes=None)`. The CLI `download` command defaults to CSV/HTML unless you specify particular suffixes or explicit `--file` values.