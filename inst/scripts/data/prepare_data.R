################################################################################
# Prepare Benchmark Data
# Scripts for preparing real and simulated benchmark datasets
################################################################################

library(dplyr)
library(Matrix)

################################################################################
# Data Structure
################################################################################

# Each benchmark dataset should be an RDS file containing a list with:
#
# $bulk           - Bulk expression matrix (genes x samples)
# $cellTypeExpr   - Reference expression matrix (genes x cell types)
# $bulkRatio      - Ground truth proportions (cell types x samples)
# $singleCellExpr - Single-cell expression matrix (genes x cells) [optional]
# $cell_type_mapping - Mapping of cell type names [optional]

################################################################################
# Example: Create Simulated Dataset
################################################################################

create_simulated_dataset <- function(n_genes = 5000,
                                      n_cell_types = 10,
                                      n_samples = 50,
                                      n_cells_per_type = 100,
                                      noise_level = 0.1,
                                      seed = 42) {
  
  set.seed(seed)
  
  message("Creating simulated dataset...")
  message(sprintf("  Genes: %d, Cell types: %d, Samples: %d", n_genes, n_cell_types, n_samples))
  
  # Generate cell-type specific gene expression profiles
  # Use negative binomial distribution to simulate count data
  base_expression <- matrix(
    rnbinom(n_genes * n_cell_types, size = 2, mu = 50),
    nrow = n_genes,
    ncol = n_cell_types
  )
  
  # Add cell-type specific marker genes
  n_markers_per_type <- 50
  for (ct in 1:n_cell_types) {
    marker_genes <- ((ct - 1) * n_markers_per_type + 1):(ct * n_markers_per_type)
    marker_genes <- marker_genes[marker_genes <= n_genes]
    base_expression[marker_genes, ct] <- base_expression[marker_genes, ct] * 5
  }
  
  rownames(base_expression) <- paste0("Gene", 1:n_genes)
  colnames(base_expression) <- paste0("CellType", 1:n_cell_types)
  
  # Generate single-cell expression (optional, for methods that require it)
  n_cells_total <- n_cells_per_type * n_cell_types
  single_cell_expr <- matrix(0, nrow = n_genes, ncol = n_cells_total)
  cell_labels <- rep(1:n_cell_types, each = n_cells_per_type)
  
  for (i in 1:n_cells_total) {
    ct <- cell_labels[i]
    single_cell_expr[, i] <- rnbinom(n_genes, size = 2, mu = base_expression[, ct])
  }
  
  rownames(single_cell_expr) <- rownames(base_expression)
  colnames(single_cell_expr) <- paste0("Cell", 1:n_cells_total)
  
  # Reference expression = average per cell type
  reference <- sapply(1:n_cell_types, function(ct) {
    cells <- which(cell_labels == ct)
    rowMeans(single_cell_expr[, cells])
  })
  colnames(reference) <- colnames(base_expression)
  
  # Generate true proportions (Dirichlet-like)
  true_proportions <- matrix(
    rgamma(n_samples * n_cell_types, shape = 1, rate = 1),
    nrow = n_cell_types,
    ncol = n_samples
  )
  true_proportions <- sweep(true_proportions, 2, colSums(true_proportions), "/")
  rownames(true_proportions) <- colnames(reference)
  colnames(true_proportions) <- paste0("Sample", 1:n_samples)
  
  # Generate bulk expression
  bulk <- reference %*% true_proportions
  
  # Add noise
  noise <- matrix(rnorm(n_genes * n_samples, sd = noise_level * mean(bulk)), 
                  nrow = n_genes, ncol = n_samples)
  bulk <- bulk + noise
  bulk[bulk < 0] <- 0
  
  colnames(bulk) <- colnames(true_proportions)
  
  # Create output list
  dataset <- list(
    bulk = bulk,
    cellTypeExpr = reference,
    bulkRatio = true_proportions,
    singleCellExpr = single_cell_expr,
    cell_type_mapping = data.frame(
      original = colnames(reference),
      clean = colnames(reference)
    )
  )
  
  message("Dataset created successfully.")
  
  return(dataset)
}

################################################################################
# Example: Prepare Real Dataset from GEO
################################################################################

prepare_geo_dataset <- function(geo_accession,
                                 reference_path,
                                 ground_truth_path,
                                 output_path) {
  
  message(sprintf("Preparing dataset: %s", geo_accession))
  
  # This is a template - actual implementation depends on data source
  # 
  # Steps:
  # 1. Download bulk RNA-seq data from GEO
  # 2. Load reference expression profiles
  # 3. Load ground truth proportions (e.g., from flow cytometry)
  # 4. Align gene names and sample IDs
  # 5. Save as standardized RDS file
  
  stop("This is a template function. Implement for your specific dataset.")
}

################################################################################
# Validate Dataset Structure
################################################################################

validate_dataset <- function(dataset) {
  
  required_fields <- c("bulk", "cellTypeExpr", "bulkRatio")
  
  # Check required fields
  for (field in required_fields) {
    if (!(field %in% names(dataset))) {
      stop(sprintf("Missing required field: %s", field))
    }
  }
  
  # Check dimensions
  bulk <- dataset$bulk
  reference <- dataset$cellTypeExpr
  ground_truth <- dataset$bulkRatio
  
  # Genes should match between bulk and reference
  common_genes <- intersect(rownames(bulk), rownames(reference))
  if (length(common_genes) < 100) {
    warning(sprintf("Only %d common genes between bulk and reference", length(common_genes)))
  }
  
  # Cell types should match between reference and ground truth
  if (!all(colnames(reference) %in% rownames(ground_truth))) {
    warning("Cell type mismatch between reference and ground truth")
  }
  
  # Samples should match between bulk and ground truth
  if (!all(colnames(bulk) %in% colnames(ground_truth))) {
    warning("Sample mismatch between bulk and ground truth")
  }
  
  # Ground truth should sum to ~1
  col_sums <- colSums(ground_truth)
  if (any(abs(col_sums - 1) > 0.01)) {
    warning("Ground truth proportions do not sum to 1 for some samples")
  }
  
  message("Dataset validation complete.")
  
  return(TRUE)
}

################################################################################
# Main Execution
################################################################################

if (sys.nframe() == 0) {
  
  # Create example simulated dataset
  sim_data <- create_simulated_dataset(
    n_genes = 5000,
    n_cell_types = 10,
    n_samples = 50
  )
  
  # Validate
  validate_dataset(sim_data)
  
  # Save
  output_dir <- "data/simulated_datasets"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(sim_data, file.path(output_dir, "example_simulated.rds"))
  message("Saved to: ", file.path(output_dir, "example_simulated.rds"))
}
