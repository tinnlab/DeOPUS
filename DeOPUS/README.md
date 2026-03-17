# DeOPUS: Deconvolution via Optimized Power-transformed Unmixing with Shrinkage

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DeOPUS is a reference-based cellular deconvolution method that employs hierarchical shrinkage transformation to robustly estimate cell-type proportions from bulk RNA sequencing data.

## Overview

Single-cell RNA sequencing enables comprehensive transcriptomic profiling at single-cell resolution, but high costs limit its widespread application. DeOPUS offers a cost-effective alternative by computationally estimating cell-type proportions from bulk RNA-seq data.

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
install.packages(c("Matrix", "dplyr", "parallel"))
```

## Quick Start

```r
library(DeOPUS)

# Load your data
# cellTypeExpr: reference expression matrix (genes x cell types)
# bulk: bulk expression matrix (genes x samples)

# Run deconvolution
results <- deconvolve(
  bulk = bulk_matrix,
  reference = reference_matrix,
  alpha = 0.01,
  n_cores = 4
)

# View estimated proportions
head(results$proportions)
```

## Input Data Format

DeOPUS requires two main inputs:

1. **Reference expression matrix** (`reference`): A genes × cell types matrix containing average expression profiles for each cell type, typically derived from scRNA-seq data.

2. **Bulk expression matrix** (`bulk`): A genes × samples matrix containing bulk RNA-seq expression data to be deconvolved.

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

The optimization minimizes the weighted loss between transformed bulk and reconstructed expression profiles using L-BFGS-B with box constraints.

## Benchmarking

We benchmarked DeOPUS against six state-of-the-art methods:
- MuSiC
- AutoGeneS
- CIBERSORT
- FARDEEP
- scaden
- AdRoit

### Running Benchmarks

```r
# Run benchmark on simulated data
source("scripts/benchmark/run_benchmark_simulated.R")

# Run benchmark on real data
source("scripts/benchmark/run_benchmark_real.R")

# Generate visualization dashboard
source("scripts/analysis/visualize_results.R")
results <- create_summary_barplot_dashboard(benchmark_data)
```

## Reproducing Paper Results

To reproduce the results from our paper:

```bash
# Clone the repository
git clone https://github.com/tinnlab/DeOPUS.git
cd DeOPUS

# Run the complete benchmark pipeline
Rscript scripts/benchmark/run_all_benchmarks.R

# Generate figures
Rscript scripts/analysis/generate_figures.R
```

### Data Availability

Benchmark datasets are available at https://doi.org/10.5281/zenodo.19050845 or can be generated using:

```r
source("scripts/data/prepare_benchmark_data.R")
```

## Output

DeOPUS returns a list containing:

- `proportions`: Matrix of estimated cell-type proportions (samples × cell types)
- `convergence`: Optimization convergence information for each sample

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 0.01 | Regularization parameter for transformation |
| `n_cores` | 1 | Number of parallel cores |
| `maxit` | 100 | Maximum optimization iterations |
| `power` | 2 | Loss function power (1=MAE, 2=MSE) |

## Citation

If you use DeOPUS in your research, please cite:

```bibtex
@article{DeOPUS2025,
  title={DeOPUS: Deconvolution via Optimized Power-transformed Unmixing with Shrinkage},
  author={},
  journal={},
  year={2025},
  doi={}
}
```

[//]: # (## License)

[//]: # ()
[//]: # (This project is licensed under the MIT License - see the [LICENSE]&#40;LICENSE&#41; file for details.)

[//]: # ()
[//]: # (## Contact)

[//]: # ()
[//]: # (For questions or issues, please open an issue on GitHub or contact [email].)
