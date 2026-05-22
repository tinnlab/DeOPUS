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

### Dependencies

```r
install.packages(c("Matrix", "parallel"))
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

Benchmark and visualization scripts ship in `inst/scripts/`. After installation,
source them via `system.file()`:

```r
# Run benchmark on real datasets
source(system.file("scripts/benchmark/run_benchmark_real.R", package = "DeOPUS"))
run_benchmark()

# Generate visualization dashboard
source(system.file("scripts/analysis/visualize_results.R", package = "DeOPUS"))
results <- create_summary_barplot_dashboard(benchmark_data)
```

Or from a clone of the repository:

```r
source("inst/scripts/benchmark/run_benchmark_real.R")
source("inst/scripts/analysis/visualize_results.R")
```

## Reproducing Paper Results

To reproduce the results from our paper:

```bash
# Clone the repository
git clone https://github.com/tinnlab/DeOPUS.git
cd DeOPUS

# Run the benchmark pipeline
Rscript inst/scripts/benchmark/run_benchmark_real.R

# Generate figures
Rscript inst/scripts/analysis/generate_figures.R
```

### Data Availability

Benchmark datasets are available at https://doi.org/10.5281/zenodo.19050845 or can
be regenerated using:

```r
source("inst/scripts/data/prepare_data.R")
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
