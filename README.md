# DeOPUS: Deconvolution via Optimized Power-transformed Unmixing with Shrinkage

[![R-CMD-check](https://github.com/tinnlab/DeOPUS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tinnlab/DeOPUS/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DeOPUS is a reference-based cellular deconvolution method that employs hierarchical
shrinkage transformation to robustly estimate cell-type proportions from bulk RNA-seq data.

## Overview

Single-cell RNA sequencing enables comprehensive transcriptomic profiling at single-cell
resolution, but high costs limit its widespread application. DeOPUS offers a cost-effective
alternative by computationally estimating cell-type proportions from bulk RNA-seq data.

**Key features:**
- Hierarchical shrinkage transformation with local and global priors
- Variance-stabilizing power transformation
- Quantile normalization to minimize outlier influence
- Robust performance across diverse tissues and cell-type complexities

## Installation

### From GitHub

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DeOPUS
devtools::install_github("tinnlab/DeOPUS")
```

### Optional dependencies (only for the benchmark / figure scripts)

The core `deconvolve()` function needs nothing beyond what `install_github`
installs. The reproducibility scripts in `inst/scripts/` have additional needs:

```r
# Benchmark runner
install.packages(c("RhpcBLASctl", "dplyr"))

# Figure generation
install.packages(c("tidyverse", "scales", "cowplot", "gridExtra"))

# Bioconductor (used by visualize_results.R)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))

# OPTIONAL — only for competing methods (MuSiC / AutoGeneS / CIBERSORT /
# FARDEEP / scaden / AdRoit). Requires Docker installed and running.
# install.packages("DeconBenchmark")
```

## Quick Start

The package ships with a built-in simulated benchmark dataset (`sampleData`):
11,852 genes × 512 samples with known ground-truth proportions for 2 cell types.

```r
library(DeOPUS)

# Load the built-in sample benchmark dataset
data(sampleData)

# Run deconvolution on a small subset for speed
set.seed(42)
idx <- sample(ncol(sampleData$bulk), 10)

results <- deconvolve(
  bulk      = sampleData$bulk[, idx],
  reference = sampleData$cellTypeExpr,
  alpha     = 0.01,
  n_cores   = 1,
  verbose   = TRUE
)

# View estimated proportions (samples x cell types)
head(results$proportions)

# Evaluate against ground truth
# IMPORTANT: bulkRatio rows and cellTypeExpr columns may be in different orders,
# so align by name before computing correlations.
ct <- colnames(results$proportions)
true_props <- t(sampleData$bulkRatio[ct, idx, drop = FALSE])  # 10 x 2 aligned

cor_values <- sapply(seq_len(nrow(results$proportions)), function(i) {
  cor(results$proportions[i, ], true_props[i, ], method = "pearson")
})
cat("Mean Pearson correlation:", round(mean(cor_values, na.rm = TRUE), 3), "\n")
```

## Input Data Format

DeOPUS requires two main inputs:

1. **Reference expression matrix** (`reference`): A genes × cell types matrix
   containing average expression profiles for each cell type, typically derived
   from scRNA-seq data.
2. **Bulk expression matrix** (`bulk`): A genes × samples matrix containing bulk
   RNA-seq expression data to be deconvolved.

Both matrices should:
- Have matching gene identifiers (rownames)
- Be in linear scale (not log-transformed)
- Contain non-negative values

## Method Details

DeOPUS applies a multi-level adaptive transformation:

1. **Local shrinkage**: Attenuates the influence of individual high-variance genes
2. **Global shrinkage**: Provides overall regularization across the expression profile
3. **Power transformation**: Stabilizes variance across the dynamic range
4. **Quantile normalization**: Ensures robust comparison between predicted and observed profiles

The optimization minimizes the weighted loss between transformed bulk and reconstructed
expression profiles using L-BFGS-B with box constraints.

## Benchmarking

We benchmarked DeOPUS against six state-of-the-art methods:
- MuSiC
- AutoGeneS
- CIBERSORT
- FARDEEP
- scaden
- AdRoit

### Running Benchmarks

Benchmark and visualization scripts ship in `inst/scripts/` and are accessible
after install via `system.file()`. The benchmark expects a directory of `.rds`
input files where each file is a list matching the schema of `data(sampleData)`
(`$bulk`, `$cellTypeExpr`, `$bulkRatio`, etc.).

To benchmark **DeOPUS** on the bundled `sampleData` end-to-end:

```r
library(DeOPUS)
source(system.file("scripts/benchmark/run_benchmark_real.R", package = "DeOPUS"))

# Prepare an input dir with at least one dataset .rds file
input_dir  <- file.path(tempdir(), "decopus_inputs")
output_dir <- file.path(tempdir(), "decopus_results")
dir.create(input_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the bundled sampleData as the benchmark input
data(sampleData)
saveRDS(sampleData, file.path(input_dir, "sampleData.rds"))

# Run benchmark — restrict to DeOPUS unless you have DeconBenchmark + Docker
run_benchmark(
  data_dir   = input_dir,
  output_dir = output_dir,
  methods    = "DeOPUS"
)

# Aggregate per-sample, per-cell-type results into a long-format dataframe
benchmark_data <- aggregate_results(output_dir = output_dir, methods = "DeOPUS")
head(benchmark_data)
```

To include MuSiC, AutoGeneS, CIBERSORT, FARDEEP, scaden, AdRoit, install the
optional `DeconBenchmark` package (requires Docker) and pass them in
`methods = c("DeOPUS", "MuSiC", ...)`.

## Reproducing Paper Results

The benchmark scripts expect a directory of `.rds` datasets (one per real
benchmark dataset, matching the schema of `data(sampleData)`) at
`data/real_datasets/`, and write results to `results/real_benchmark/` and
figures to `figures/`.

```bash
git clone https://github.com/tinnlab/DeOPUS.git
cd DeOPUS

# Place benchmark datasets here (one .rds per dataset)
mkdir -p data/real_datasets
# cp /path/to/your/*.rds data/real_datasets/

# Run benchmark pipeline (defaults: data/real_datasets/ → results/real_benchmark/)
Rscript inst/scripts/benchmark/run_benchmark_real.R

# Generate figures (written to figures/)
Rscript inst/scripts/analysis/generate_figures.R
```

### Data Availability

Benchmark datasets are available at
https://doi.org/10.5281/zenodo.19050845. The simulation generator and a
template for preparing real GEO datasets are in
`inst/scripts/data/prepare_data.R` — source the file, then call
`create_simulated_dataset()` or adapt `prepare_geo_dataset()` to your data:

```r
source(system.file("scripts/data/prepare_data.R", package = "DeOPUS"))
# sampleData_like <- create_simulated_dataset(n_samples = 100, n_cell_types = 5)
```

## Output

`deconvolve()` returns a list containing:

- `proportions`: Matrix of estimated cell-type proportions (samples × cell types). Each row sums to 1.
- `convergence`: Named list of per-sample optimization convergence information (one entry per sample, holding `convergence`, `value`, and optionally an `error` message).

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha`   | 0.01    | Regularization parameter for transformation |
| `power`   | 2       | Loss function power (1 = MAE, 2 = MSE) |
| `n_cores` | 1       | Number of parallel cores |
| `maxit`   | 100     | Maximum optimization iterations |
| `verbose` | FALSE   | Print progress messages |

## Citation

If you use DeOPUS in your research, please cite:

```bibtex
@article{DeOPUS2025,
  title  = {DeOPUS: Deconvolution via Optimized Power-transformed Unmixing with Shrinkage},
  author = {Ha Nguyen},
  year   = {2025}
}
```

## License

MIT © Ha Nguyen — see [LICENSE.md](LICENSE.md) for full text.
