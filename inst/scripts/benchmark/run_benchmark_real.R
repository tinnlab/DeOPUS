################################################################################
# Benchmark Script for Real Datasets
# Compares DeOPUS against competing deconvolution methods
################################################################################

# Set thread limits to avoid oversubscription
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

# Load required packages
library(Matrix)
library(dplyr)
library(parallel)
library(DeOPUS)

# Optional: Load DeconBenchmark for competing methods
# library(DeconBenchmark)
# library(babelwhale)

################################################################################
# Configuration
################################################################################

# Define methods to benchmark
METHODS <- c("DeOPUS", "MuSiC", "AutoGeneS", "CIBERSORT", "FARDEEP", "scaden", "AdRoit")

# Directories (modify as needed)
DATA_DIR <- "data/real_datasets"
OUTPUT_DIR <- "results/real_benchmark"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Number of cores for parallel processing
N_CORES <- parallel::detectCores() - 1

################################################################################
# DeOPUS Deconvolution Function
################################################################################

run_decopus <- function(data, n_cores = 1) {
  
  cellTypeExpr <- as.matrix(data$cellTypeExpr)
  bulk <- as.matrix(data$bulk)
  
  # Run DeOPUS
  results <- deconvolve(
    bulk = bulk,
    reference = cellTypeExpr,
    alpha = 0.01,
    power = 2,
    n_cores = n_cores,
    maxit = 100
  )
  
  # Transpose to match expected format (cell types x samples)
  res <- t(results$proportions)
  
  return(res)
}

################################################################################
# Run Benchmark for Competing Methods (requires DeconBenchmark)
################################################################################

run_competing_method <- function(data, method, docker_args = NULL) {
  
  # Prepare data for DeconBenchmark
  benchmark_data <- list(
    bulk = as.matrix(data$bulk),
    cellTypeExpr = as.matrix(data$cellTypeExpr),
    singleCellExpr = as.matrix(data$singleCellExpr)
  )
  
  # Default Docker arguments
  if (is.null(docker_args)) {
    docker_args <- c(
      '--cpus=2.0',
      '-m=16G',
      '--memory-reservation=12G',
      '-e OMP_NUM_THREADS=1',
      '-e OPENBLAS_NUM_THREADS=1',
      '-e MKL_NUM_THREADS=1'
    )
  }
  
  # Run deconvolution using DeconBenchmark
  res <- DeconBenchmark::runDeconvolution(
    methods = method,
    verbose = TRUE,
    dockerArgs = docker_args,
    timeout = 12 * 3600,
    bulk = benchmark_data$bulk,
    cellTypeExpr = benchmark_data$cellTypeExpr,
    singleCellExpr = benchmark_data$singleCellExpr
  )
  
  return(res[[method]]$P)
}

################################################################################
# Main Benchmark Loop
################################################################################

run_benchmark <- function(data_dir = DATA_DIR,
                          output_dir = OUTPUT_DIR,
                          methods = METHODS,
                          n_cores = N_CORES) {
  
  # List all data files
  data_files <- list.files(data_dir, pattern = "\\.rds$", full.names = FALSE)
  
  message(sprintf("Found %d datasets to benchmark", length(data_files)))
  message(sprintf("Methods: %s", paste(methods, collapse = ", ")))
  
  # Process each dataset
  for (file in data_files) {
    message(sprintf("\n=== Processing: %s ===", file))
    
    # Load data
    data_path <- file.path(data_dir, file)
    data <- readRDS(data_path)
    
    if (is.null(data)) {
      warning(sprintf("Failed to load: %s", file))
      next
    }
    
    # Extract ground truth
    ground_truth <- data$bulkRatio
    cell_type_mapping <- data$cell_type_mapping
    
    # Run each method
    for (method in methods) {
      message(sprintf("  Running %s...", method))
      
      # Create output path
      method_dir <- file.path(output_dir, method)
      if (!dir.exists(method_dir)) {
        dir.create(method_dir, recursive = TRUE)
      }
      
      res_file <- file.path(method_dir, file)
      
      # Skip if already processed
      if (file.exists(res_file)) {
        message(sprintf("    Skipping (already exists)"))
        next
      }
      
      # Run deconvolution
      tryCatch({
        if (method == "DeOPUS") {
          res <- run_decopus(data, n_cores = n_cores)
        } else {
          res <- run_competing_method(data, method)
        }
        
        # Align samples
        if (!is.null(res)) {
          samples_common <- intersect(colnames(ground_truth), colnames(res))
          if (length(samples_common) == 0) {
            res <- t(res)
            samples_common <- intersect(colnames(ground_truth), colnames(res))
          }
          
          # Save results
          saveRDS(list(
            res = res,
            grountTruth = ground_truth,
            cell_type_mapping = cell_type_mapping
          ), res_file)
          
          message(sprintf("    Saved to: %s", res_file))
        }
        
      }, error = function(e) {
        warning(sprintf("    Error: %s", conditionMessage(e)))
      })
    }
  }
  
  message("\n=== Benchmark Complete ===")
}

################################################################################
# Calculate Metrics
################################################################################

calculate_metrics <- function(res, ground_truth) {
  
  # Ensure matching dimensions
  common_cells <- intersect(rownames(res), rownames(ground_truth))
  common_samples <- intersect(colnames(res), colnames(ground_truth))
  
  res <- res[common_cells, common_samples, drop = FALSE]
  ground_truth <- ground_truth[common_cells, common_samples, drop = FALSE]
  
  # Calculate metrics per sample (column)
  metrics <- list(
    cor_pearson = sapply(seq_len(ncol(res)), function(i) {
      cor(res[, i], ground_truth[, i], method = "pearson", use = "complete.obs")
    }) %>% mean(na.rm = TRUE),
    
    cor_spearman = sapply(seq_len(ncol(res)), function(i) {
      cor(res[, i], ground_truth[, i], method = "spearman", use = "complete.obs")
    }) %>% mean(na.rm = TRUE),
    
    mse = sapply(seq_len(ncol(res)), function(i) {
      mean((res[, i] - ground_truth[, i])^2)
    }) %>% mean(na.rm = TRUE),
    
    mae = sapply(seq_len(ncol(res)), function(i) {
      mean(abs(res[, i] - ground_truth[, i]))
    }) %>% mean(na.rm = TRUE)
  )
  
  return(metrics)
}

################################################################################
# Aggregate Results
################################################################################

aggregate_results <- function(output_dir = OUTPUT_DIR, methods = METHODS) {
  
  all_results <- list()
  
  for (method in methods) {
    method_dir <- file.path(output_dir, method)
    
    if (!dir.exists(method_dir)) {
      next
    }
    
    result_files <- list.files(method_dir, pattern = "\\.rds$", full.names = TRUE)
    
    for (res_file in result_files) {
      dataset <- basename(res_file)
      
      tryCatch({
        data <- readRDS(res_file)
        
        res <- data$res
        ground_truth <- data$grountTruth
        
        # Align cell types
        common_cells <- intersect(rownames(res), rownames(ground_truth))
        common_samples <- intersect(colnames(res), colnames(ground_truth))
        
        if (length(common_cells) == 0 || length(common_samples) == 0) {
          next
        }
        
        res <- res[common_cells, common_samples, drop = FALSE]
        ground_truth <- ground_truth[common_cells, common_samples, drop = FALSE]
        
        # Create long-format data for each sample and cell type
        for (sample in common_samples) {
          for (cell_type in common_cells) {
            all_results[[length(all_results) + 1]] <- data.frame(
              method = method,
              dataset = dataset,
              sample = sample,
              cell_type = cell_type,
              ground_truth = ground_truth[cell_type, sample],
              predicted_value = res[cell_type, sample],
              stringsAsFactors = FALSE
            )
          }
        }
        
      }, error = function(e) {
        warning(sprintf("Error processing %s/%s: %s", method, dataset, conditionMessage(e)))
      })
    }
  }
  
  # Combine all results
  results_df <- do.call(rbind, all_results)
  
  return(results_df)
}

################################################################################
# Main Execution
################################################################################

if (sys.nframe() == 0) {
  # Run benchmark
  run_benchmark()
  
  # Aggregate results
  all_data <- aggregate_results()
  
  # Save aggregated data
  saveRDS(all_data, file.path(OUTPUT_DIR, "benchmark_results_aggregated.rds"))
  
  message("Results saved to: ", file.path(OUTPUT_DIR, "benchmark_results_aggregated.rds"))
}
