<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
</head>
<body>

<h1>CIBERSORT (LM22) Runner — Python + rpy2</h1>

<p>Run LM22-based immune cell deconvolution (CIBERSORT) from Python via <strong>rpy2</strong>, with an automated <strong>QA pass</strong> that validates and labels samples (“Excellent / Moderate / Poor”).</p>

<blockquote>
<p>This runner executes an R pipeline under the hood, handles minimal preprocessing (counts→CPM via edgeR when needed), writes standard plots/CSVs, and performs a numeric-safety QA step.</p>
</blockquote>

<hr />

<h2>Features</h2>

<ul>
<li><strong>Input auto-detection</strong>: Treats inputs as TPM/RPKM if sample sums ≈ 1e6; otherwise computes <strong>TMM-normalized CPM</strong> (edgeR).</li>
<li><strong>Metadata-aware</strong>: Optional metadata alignment by <code>sample_id</code>.</li>
<li><strong>Artifacts</strong>:
  <ul>
    <li><code>CIBERSORT_results.csv</code> (fractions + metrics)</li>
    <li>Stacked barplots per sample (paged)</li>
    <li>Heatmap (samples × cell types)</li>
    <li>P-value histogram; Corr vs RMSE scatter</li>
    <li>LM22 overlap reports</li>
    <li><code>CIBERSORT_Quality_Assessment.csv</code> with categorical labels</li>
    <li><code>RUN_SUMMARY.txt</code> with session info</li>
  </ul>
</li>
<li><strong>Fail-fast</strong> for missing inputs and non-numeric columns.</li>
</ul>

<hr />

<h2>Requirements</h2>

<ul>
<li><strong>Python</strong> 3.9–3.12</li>
<li><strong>R</strong> ≥ 4.0 installed and on <code>PATH</code></li>
<li><strong>Python packages</strong>: <code>rpy2</code></li>
<li><strong>R packages</strong> (installed automatically if missing):
  <ul>
    <li>CRAN: <code>devtools</code>, <code>readr</code>, <code>readxl</code>, <code>dplyr</code>, <code>tibble</code>, <code>stringr</code>, <code>tools</code>, <code>ggplot2</code>, <code>pheatmap</code>, <code>reshape2</code>, <code>tidyr</code>, <code>ggrepel</code>, <code>scales</code>, <code>data.table</code></li>
    <li>Bioc: <code>edgeR</code></li>
    <li>GitHub (if needed): <code>Moonerss/CIBERSORT</code> (installed via <code>devtools::install_github</code>)</li>
  </ul>
</li>
</ul>

<p>The runner also looks for a local folder <code>CIBERSORT-main</code> and will <code>install_local</code> if present.</p>

<hr />

<h2>Install</h2>

<h3>Conda (recommended)</h3>
<pre><code>conda create -n cibersort python=3.11 -y
conda activate cibersort
# System R (Conda build)
conda install -c conda-forge r-base -y
# Python deps
pip install rpy2
</code></pre>

<p>If using a system R outside Conda, ensure <code>R</code> is discoverable: <code>R --version</code> works in the same shell.</p>

<hr />

<h2>CLI</h2>
<pre><code>python Cibersort.py --counts &lt;file&gt; --lm22 &lt;file&gt; --out &lt;dir&gt; [options]
</code></pre>

<h3>Required</h3>
<ul>
<li><code>--counts</code> : TSV/CSV/XLS(X). First column = gene symbol (header names are auto-detected).</li>
<li><code>--lm22</code>   : Path to LM22 signature matrix (txt/tsv).</li>
<li><code>--out</code>    : Output directory.</li>
</ul>

<h3>Optional</h3>
<ul>
<li><code>--meta</code>         : Metadata file with columns <code>sample_id, condition</code>.</li>
<li><code>--perm</code>         : CIBERSORT permutations (default: 100).</li>
<li><code>--qn</code>           : Enable quantile normalization (<code>true/false</code>; default false; keep <code>false</code> for RNA-seq/TPM).</li>
<li><code>--chunk-size</code>   : Samples per page for stacked bar plot (default: 60).</li>
<li><code>--install</code>      : Pre-install base R packages before running.</li>
<li><code>--res-path</code>     : Explicit path to <code>CIBERSORT_results.csv</code> for the QA step (rarely needed).</li>
<li><code>--log-level</code>    : <code>DEBUG|INFO|WARNING|ERROR</code> (default: <code>INFO</code>).</li>
</ul>

<h3>Exit codes</h3>
<ul>
<li><code>0</code> success</li>
<li><code>1</code> pipeline exception</li>
<li><code>2</code> missing input files</li>
</ul>

<hr />

<h2>Usage Examples</h2>

<h3>Bash / sh (Linux/macOS)</h3>
<pre><code>python Cibersort.py \
  --counts 'Bulk-data/counts.csvv' \
  --meta   'Bulk-data/meta.csv' \
  --lm22   'inst/extdata/LM22.txt' \
  --out    'CIBERSORT_outputs-v' \
  --perm 100 --qn false --chunk-size 40 --install
</code></pre>

<p>Paths with spaces/parentheses <strong>must be quoted</strong> (or escaped) on all shells.</p>

<h3>Windows PowerShell</h3>
<pre><code>python Cibersort.py `
  --counts "Bulk-data/counts.csvv" `
  --meta   "Bulk-data/meta.csv" `
  --lm22   "inst/extdata/LM22.txt" `
  --out    "CIBERSORT_outputs-v" `
  --perm 100 --qn false --chunk-size 40 --install
</code></pre>

<h3>Windows CMD</h3>
<pre><code>python Cibersort.py ^
  --counts "Bulk-data/counts.csvv" ^
  --meta   "Bulk-data/meta.csv" ^
  --lm22   "inst/extdata/LM22.txt" ^
  --out    "CIBERSORT_outputs-v" ^
  --perm 100 --qn false --chunk-size 40 --install
</code></pre>

<hr />

<h2>Inputs &amp; Format</h2>

<h3>Counts / Expression</h3>
<ul>
<li>Row = gene, Column = sample.</li>
<li>First column = gene symbol. Accepted headers include: <code>Gene</code>, <code>gene</code>, <code>GeneSymbol</code>, <code>Symbol</code>, <code>ID</code>, etc.</li>
<li>The runner:
  <ul>
    <li>Keeps genes with expression <code>&gt; 1</code> in ≥10% of samples (≥1 sample minimum).</li>
    <li>Auto-detects scale:
      <ul>
        <li>Median column sum ~ 1e6 → treat as TPM/RPKM (<code>QN=FALSE</code> recommended).</li>
        <li>Otherwise → edgeR TMM → CPM.</li>
      </ul>
    </li>
  </ul>
</li>
</ul>

<h3>Metadata (optional)</h3>
<ul>
<li>Columns:  
  <ul>
    <li><code>sample_id</code> (must match counts’ column names)</li>
    <li><code>condition</code> (free text; used only in QC table join, currently no group plots)</li>
  </ul>
</li>
</ul>

<h3>LM22</h3>
<ul>
<li>Text/TSV with first column as gene and 22 cell types as columns.</li>
<li>We upper-case gene symbols internally for overlap accounting.</li>
</ul>

<hr />

<h2>Outputs (key files)</h2>

<ul>
<li><code>CIBERSORT_results.csv</code> — Deconvolution results (rows=samples). Common trailing columns:
  <ul>
    <li><code>P-value</code>, <code>Correlation</code>, <code>RMSE</code>, <code>Absolute score</code> (if provided by CIBERSORT impl)</li>
  </ul>
</li>
<li><code>01_stacked_bar_fractions*.png</code> — Horizontal stacked bars of cell fractions (paged by <code>--chunk-size</code>)</li>
<li><code>02_heatmap_all_cell_types.png</code> — Heatmap (samples × LM22 types)</li>
<li><code>03_pvalue_histogram.png</code> — Per-sample P-value distribution (if present)</li>
<li><code>04_scatter_correlation_vs_RMSE.png</code> — Fit scatter (if both fields present)</li>
<li><code>LM22_overlap_gene_values_by_sample.csv</code> — Per-sample expression for LM22-overlap genes</li>
<li><code>LM22_overlap_gene_summary.csv</code> — Mean expression + cell types per gene</li>
<li><code>LM22_overlap_report.txt</code> — Coverage summary</li>
<li><code>CIBERSORT_Quality_Assessment.csv</code> — QA summary per sample:
  <ul>
    <li>Fraction sum check (0.85–1.15), discretized P/Corr/RMSE, <code>Quality_Category</code></li>
  </ul>
</li>
<li><code>RUN_SUMMARY.txt</code> — Inputs, counts, session info, medians</li>
</ul>

<hr />

<h2>Logging</h2>

<ul>
<li>Python logs to stdout (<code>--log-level DEBUG</code> for verbose).</li>
<li>R messages are surfaced via rpy2.</li>
<li>Failures during R package install or CIBERSORT run are bubbled up and return exit code <code>1</code>.</li>
</ul>

<hr />

# Python deps
RUN pip install --no-cache-dir rpy2


# Copy code
WORKDIR /app
COPY Cibersort.py .


# Default command (print help)
CMD ["python", "Cibersort.py", "--help"]
</code></pre>

<p>Build &amp; run:</p>
<pre><code>docker build -t cibersort-runner .
docker run --rm -v "$PWD":/work -w /work cibersort-runner \
  python /app/Cibersort.py --help
</code></pre>

<hr />
>

<hr />

<h2>Maintenance Checklist</h2>

<ul>
<li>Validate paths quoting on all OSes.</li>
<li>Confirm R package installs on clean hosts (<code>--install</code> first run).</li>
<li>Keep LM22 reference versioned in <code>inst/extdata/</code>.</li>
<li>CI smoke-test: run with a tiny synthetic dataset (see below) to ensure plots/CSVs render.</li>
</ul>

<h3>Tiny Synthetic Smoke Test (optional)</h3>
<pre><code># create_fake.py
import numpy as np, pandas as pd
rng = np.random.default_rng(1)
genes = [f"GENE{i}" for i in range(200)]
samples = [f"S{i}" for i in range(12)]
df = pd.DataFrame(rng.poisson(5, size=(len(genes), len(samples))), index=genes, columns=samples)
df.insert(0, "GeneSymbol", genes)
df.to_csv("fake_counts.csv", index=False)
print("fake_counts.csv written")
</code></pre>
<p>Then:</p>
<pre><code>python create_fake.py
python Cibersort.py --counts fake_counts.csv --lm22 inst/extdata/LM22.txt --out out_fake --perm 10 --qn false --install
</code></pre>

<hr />

<h2>FAQ</h2>

<ul>
<li><strong>Why is QN default false?</strong><br />
  For RNA-seq/TPM, quantile normalization is generally not recommended; we leave it opt-in.</li>
<li><strong>Where do QA thresholds come from?</strong><br />
  They’re pragmatic defaults: P≤0.05, Corr≥0.5, RMSE≤0.7 and fraction sum in [0.85, 1.15]. Tune in the R QA block if needed.</li>
<li><strong>Can I plug other signatures?</strong><br />
  Yes. Use any CIBERSORT-compatible signature matrix via <code>--lm22 &lt;path&gt;</code>; plots/QA still work.</li>
</ul>

<hr />

<h3>One-liner to remember (Linux/macOS)</h3>
<pre><code>python Cibersort.py --counts 'Bulk-data/counts.csvv' \
  --meta 'Bulk-data/meta.csv' \
  --lm22 'inst/extdata/LM22.txt' \
  --out 'CIBERSORT_outputs-v' --perm 100 --qn false --chunk-size 40 --install
</code></pre>

</body>
</html>
