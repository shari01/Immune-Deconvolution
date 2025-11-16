# CIBERSORT (LM22) Runner — Python + rpy2

Run LM22-based immune cell deconvolution (CIBERSORT) from Python via **rpy2**, with an automated **QA pass** that validates and labels samples (“Excellent / Moderate / Poor”).

> This runner executes an R pipeline under the hood, handles minimal preprocessing (counts→CPM via edgeR when needed), writes standard plots/CSVs, and performs a numeric-safety QA step.

---

## Features

- **Input auto-detection**: Treats inputs as TPM/RPKM if sample sums ≈ 1e6; otherwise computes **TMM-normalized CPM** (edgeR).
- **Metadata-aware**: Optional metadata alignment by `sample_id`.
- **Artifacts**:
  - `CIBERSORT_results.csv` (fractions + metrics)
  - Stacked barplots per sample (paged)
  - Heatmap (samples × cell types)
  - P-value histogram; Corr vs RMSE scatter
  - LM22 overlap reports
  - `CIBERSORT_Quality_Assessment.csv` with categorical labels
  - `RUN_SUMMARY.txt` with session info
- **Fail-fast** for missing inputs and non-numeric columns.

---

## Requirements

- **Python** 3.9–3.12
- **R** ≥ 4.0 installed and on `PATH`
- **Python packages**: `rpy2`
- **R packages** (installed automatically if missing):
  - CRAN: `devtools`, `readr`, `readxl`, `dplyr`, `tibble`, `stringr`, `tools`, `ggplot2`, `pheatmap`, `reshape2`, `tidyr`, `ggrepel`, `scales`, `data.table`
  - Bioc: `edgeR`
  - GitHub (if needed): `Moonerss/CIBERSORT` (installed via `devtools::install_github`)

> The runner also looks for a local folder `CIBERSORT-main` and will `install_local` if present.

---

## Install

### Conda (recommended)
```bash
conda create -n cibersort python=3.11 -y
conda activate cibersort
# System R (Conda build)
conda install -c conda-forge r-base -y
# Python deps
pip install rpy2
```

> If using a system R outside Conda, ensure `R` is discoverable: `R --version` works in the same shell.

---

## CLI

```
python Cibersort.py --counts <file> --lm22 <file> --out <dir> [options]
```

**Required**
- `--counts` : TSV/CSV/XLS(X). First column = gene symbol (header names are auto-detected).
- `--lm22`   : Path to LM22 signature matrix (txt/tsv).
- `--out`    : Output directory.

**Optional**
- `--meta`         : Metadata file with columns `sample_id, condition`.
- `--perm`         : CIBERSORT permutations (default: 100).
- `--qn`           : Enable quantile normalization (`true/false`; default false; keep `false` for RNA-seq/TPM).
- `--chunk-size`   : Samples per page for stacked bar plot (default: 60).
- `--install`      : Pre-install base R packages before running.
- `--res-path`     : Explicit path to `CIBERSORT_results.csv` for the QA step (rarely needed).
- `--log-level`    : `DEBUG|INFO|WARNING|ERROR` (default: `INFO`).

**Exit codes**
- `0` success
- `1` pipeline exception
- `2` missing input files

---

## Usage Examples

### Bash / zsh (Linux/macOS)
```bash
python Cibersort.py \
  --counts 'Bulk-data/TW-20250618_counts (4) (2).csv' \
  --meta   'Bulk-data/TW-20250618 _metadata (1).csv' \
  --lm22   'inst/extdata/LM22.txt' \
  --out    'CIBERSORT_outputs-v' \
  --perm 100 --qn false --chunk-size 40 --install
```

> Paths with spaces/parentheses **must be quoted** (or escaped) on all shells.

### Windows PowerShell
```powershell
python Cibersort.py `
  --counts "Bulk-data/TW-20250618_counts (4) (2).csv" `
  --meta   "Bulk-data/TW-20250618 _metadata (1).csv" `
  --lm22   "inst/extdata/LM22.txt" `
  --out    "CIBERSORT_outputs-v" `
  --perm 100 --qn false --chunk-size 40 --install
```

### Windows CMD
```cmd
python Cibersort.py ^
  --counts "Bulk-data/TW-20250618_counts (4) (2).csv" ^
  --meta   "Bulk-data/TW-20250618 _metadata (1).csv" ^
  --lm22   "inst/extdata/LM22.txt" ^
  --out    "CIBERSORT_outputs-v" ^
  --perm 100 --qn false --chunk-size 40 --install
```

---

## Inputs & Format

### Counts / Expression
- Row = gene, Column = sample.
- First column = gene symbol. Accepted headers include: `Gene`, `gene`, `GeneSymbol`, `Symbol`, `ID`, etc.
- The runner:
  - Keeps genes with expression `> 1` in ≥10% of samples (≥1 sample minimum).
  - Auto-detects scale:
    - Median column sum ~ 1e6 → treat as TPM/RPKM (`QN=FALSE` recommended).
    - Otherwise → edgeR TMM → CPM.

### Metadata (optional)
- Columns:  
  `sample_id` (must match counts’ column names)  
  `condition` (free text; used only in QC table join, currently no group plots)
  
### LM22
- Text/TSV with first column as gene and 22 cell types as columns.
- We upper-case gene symbols internally for overlap accounting.

---

## Outputs (key files)

- `CIBERSORT_results.csv` — Deconvolution results (rows=samples). Common trailing columns:
  - `P-value`, `Correlation`, `RMSE`, `Absolute score` (if provided by CIBERSORT impl)
- `01_stacked_bar_fractions*.png` — Horizontal stacked bars of cell fractions (paged by `--chunk-size`)
- `02_heatmap_all_cell_types.png` — Heatmap (samples × LM22 types)
- `03_pvalue_histogram.png` — Per-sample P-value distribution (if present)
- `04_scatter_correlation_vs_RMSE.png` — Fit scatter (if both fields present)
- `LM22_overlap_gene_values_by_sample.csv` — Per-sample expression for LM22-overlap genes
- `LM22_overlap_gene_summary.csv` — Mean expression + cell types per gene
- `LM22_overlap_report.txt` — Coverage summary
- `CIBERSORT_Quality_Assessment.csv` — QA summary per sample:
  - Fraction sum check (0.85–1.15), discretized P/Corr/RMSE, `Quality_Category`
- `RUN_SUMMARY.txt` — Inputs, counts, session info, medians

---

## Logging

- Python logs to stdout (`--log-level DEBUG` for verbose).
- R messages are surfaced via rpy2.
- Failures during R package install or CIBERSORT run are bubbled up and return exit code `1`.

---

## Troubleshooting

1. **`rpy2` cannot find R**
   - Ensure `R` on PATH: `R --version`.
   - On Conda: `conda install -c conda-forge r-base`.
   - On macOS: install R from CRAN pkg; relaunch terminal.

2. **Spaces/Parentheses in paths**
   - Quote paths exactly as they exist.  
     Example: `"Bulk-data/TW-20250618 _metadata (1).csv"` (note the space before `_metadata`).

3. **`Missing counts/lm22`**
   - The runner checks existence; use absolute paths or correct working directory.

4. **“These columns are not numeric after coercion” (QA step)**
   - Ensure cell-type columns in `CIBERSORT_results.csv` are numeric.
   - The QA tries to strip `%`, commas, spaces; if columns were renamed or mangled by Excel, keep headers simple.

5. **No overlapping `sample_id`**
   - Metadata will be ignored if `sample_id` doesn’t match counts’ sample columns.

6. **Poor QA everywhere**
   - Increase `--perm` (e.g., 1000).
   - Verify input scale (TPM vs raw). For raw counts, do **not** set `--qn true`.

---

## Reproducible Environment (Docker, optional)

```dockerfile
FROM mambaorg/micromamba:1.5.8
ARG DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-lc"]

# Create env with Python + R
RUN micromamba create -n cibersort -y -c conda-forge python=3.11 r-base r-essentials \
 && micromamba clean --all -y
SHELL ["micromamba", "run", "-n", "cibersort", "/bin/bash", "-lc"]

# Python deps
RUN pip install --no-cache-dir rpy2

# Copy code
WORKDIR /app
COPY Cibersort.py .

# Default command (print help)
CMD ["python", "Cibersort.py", "--help"]
```

Build & run:
```bash
docker build -t cibersort-runner .
docker run --rm -v "$PWD":/work -w /work cibersort-runner \
  python /app/Cibersort.py --help
```

---

## Development Notes

- **Structure**
  - Python orchestrates environment → calls two R blocks: **main pipeline** then **QA**.
  - Inputs are passed to R via `Sys.setenv` (keys prefixed `PY_…`).
- **Extending**
  - Add new R blocks and set their env vars in `run_pipeline`.
  - Keep new outputs in a dedicated subfolder to avoid name collisions.
  - Prefer writing long-form/tidy CSVs alongside wide matrices.

---

## Maintenance Checklist

- [ ] Validate paths quoting on all OSes.
- [ ] Confirm R package installs on clean hosts (`--install` first run).
- [ ] Keep LM22 reference versioned in `inst/extdata/`.
- [ ] CI smoke-test: run with a tiny synthetic dataset (see below) to ensure plots/CSVs render.

### Tiny Synthetic Smoke Test (optional)
```python
# create_fake.py
import numpy as np, pandas as pd
rng = np.random.default_rng(1)
genes = [f"GENE{i}" for i in range(200)]
samples = [f"S{i}" for i in range(12)]
df = pd.DataFrame(rng.poisson(5, size=(len(genes), len(samples))), index=genes, columns=samples)
df.insert(0, "GeneSymbol", genes)
df.to_csv("fake_counts.csv", index=False)
print("fake_counts.csv written")
```
Then:
```bash
python create_fake.py
python Cibersort.py --counts fake_counts.csv --lm22 inst/extdata/LM22.txt --out out_fake --perm 10 --qn false --install
```

---

## FAQ

- **Why is QN default false?**  
  For RNA-seq/TPM, quantile normalization is generally not recommended; we leave it opt-in.

- **Where do QA thresholds come from?**  
  They’re pragmatic defaults: P≤0.05, Corr≥0.5, RMSE≤0.7 and fraction sum in [0.85, 1.15]. Tune in the R QA block if needed.

- **Can I plug other signatures?**  
  Yes. Use any CIBERSORT-compatible signature matrix via `--lm22 <path>`; plots/QA still work.

---

## License & Ownership

- Internal research tooling. Add an explicit license header if distributing outside your org.

---

### One-liner to remember (Linux/macOS)
```bash
python Cibersort.py --counts 'Bulk-data/TW-20250618_counts (4) (2).csv' \
  --meta 'Bulk-data/TW-20250618 _metadata (1).csv' \
  --lm22 'inst/extdata/LM22.txt' \
  --out 'CIBERSORT_outputs-v' --perm 100 --qn false --chunk-size 40 --install
```
