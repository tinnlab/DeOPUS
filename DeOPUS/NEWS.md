# DeOPUS 1.0.0

## Initial Release

* First public release of DeOPUS (Deconvolution via Optimized Power-transformed Unmixing with Shrinkage)
* Core deconvolution algorithm with hierarchical shrinkage transformation
* Parallel processing support via `parallel::mclapply`
* Comprehensive benchmarking scripts for real and simulated datasets
* Visualization dashboard for benchmark results
* Full documentation and vignettes

## Key Features

* **Hierarchical shrinkage transformation**: Combines local and global shrinkage priors for robust estimation
* **Variance-stabilizing power transformation**: Handles high dynamic range in gene expression data
* **Quantile normalization**: Minimizes influence of outliers and high-variance genes
* **L-BFGS-B optimization**: Efficient constrained optimization with box constraints

## Benchmark Results

* Evaluated on 122 simulated tissues with up to 40 cell types
* Validated on 17 real bulk RNA-seq datasets with experimental ground truth
* Achieved highest Pearson correlation (mean r = 0.82) among all methods tested
* Only method to achieve positive correlations across all real datasets
