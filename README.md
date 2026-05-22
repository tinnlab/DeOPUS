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
- Parallel processing support via `parallel::mclapply`

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DeOPUS from GitHub
devtools::install_github("tinnlab/DeOPUS")
```

## Quick Start: Simulated Dataset

The package ships with a built-in simulated benchmark dataset (`simulated`) containing
11,852 genes, 512 samples, and 2 cell types with known ground-truth proportions.

```r
library(DeOPUS)

# Load the built-in simulated benchmark dataset
data(simulated)

# Run deconvolution on a small subset for speed
set.seed(42)
idx <- sample(ncol(simulated$bulk), 10)

results <- deconvolve(
  bulk      = simulated$bulk[, idx],
  reference = simulated$cellTypeExpr,
  n_cores   = 1,
  verbose   = TRUE
)

# View estimated proportions (samples x cell types)
head(results$proportions)

# Evaluate against ground truth
true_props <- t(simulated$bulkRatio[, idx])  # samples x cell types
cor_values <- sapply(seq_len(nrow(results$proportions)), function(i) {
  cor(results$proportions[i, ], true_props[i, ], method = "pearson")
})
cat("Mean Pearson correlation:", round(mean(cor_values, na.rm = TRUE), 3), "\n")
```

## Real Dataset: GSE77343

A real bulk RNA-seq dataset (GSE77343, 17,275 genes × 197 samples, 15 cell types) is
bundled in `inst/extdata/`. Load it with `system.file()`:

```r
library(DeOPUS)

# Load the real GEO dataset
GSE77343 <- readRDS(system.file("extdata", "GSE77343.rds", package = "DeOPUS"))

# Inspect
dim(GSE77343$bulk)          # 17275 x 197
dim(GSE77343$cellTypeExpr)  # 17275 x 15

# Run deconvolution
results <- deconvolve(
  bulk      = GSE77343$bulk,
  reference = GSE77343$cellTypeExpr,
  n_cores   = 4,
  verbose   = TRUE
)

# View estimated proportions
head(results$proportions)

# Compare with available ground truth (5 measurable cell types)
common_types <- intersect(rownames(GSE77343$bulkRatio), colnames(results$proportions))
cor_per_type <- sapply(common_types, function(ct) {
  cor(results$proportions[, ct], GSE77343$bulkRatio[ct, ], method = "pearson")
})
print(round(cor_per_type, 3))
```

## Input Data Format

| Argument    | Format                    | Description |
|-------------|---------------------------|-------------|
| `bulk`      | genes × samples matrix    | Bulk RNA-seq expression (linear scale, non-negative) |
| `reference` | genes × cell types matrix | Reference profiles derived from scRNA-seq |

Both matrices must share gene identifiers as row names.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha`   | 0.01    | Regularization parameter for transformation |
| `power`   | 2       | Loss function power (1 = MAE, 2 = MSE) |
| `n_cores` | 1       | Number of parallel cores |
| `maxit`   | 100     | Maximum optimization iterations |
| `verbose` | FALSE   | Print progress messages |

## Output

`deconvolve()` returns a list:

- `proportions`: Numeric matrix (samples × cell types) of estimated proportions that sum to 1.
- `convergence`: List of per-sample optimization convergence information.

## Method Details

DeOPUS applies a multi-level adaptive transformation:

1. **Local shrinkage**: Attenuates the influence of individual high-variance genes
2. **Global shrinkage**: Provides overall regularization across the expression profile
3. **Power transformation**: Stabilizes variance across the dynamic range
4. **Quantile normalization**: Ensures robust comparison between predicted and observed profiles

The optimization minimizes the weighted loss using L-BFGS-B with box constraints.

## Benchmarking Scripts

Scripts for reproducing benchmark results are installed with the package:

```r
system.file("scripts", package = "DeOPUS")
```

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

MIT © Ha Nguyen
