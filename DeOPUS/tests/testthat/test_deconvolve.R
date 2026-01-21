library(testthat)
library(DeOPUS)

test_that("deconvolve works with simple input", {
  
  # Create simple test data
  set.seed(42)
  n_genes <- 100
  n_cell_types <- 3
  n_samples <- 5
  
  # Reference matrix
  reference <- matrix(
    rpois(n_genes * n_cell_types, lambda = 100),
    nrow = n_genes,
    ncol = n_cell_types
  )
  rownames(reference) <- paste0("Gene", 1:n_genes)
  colnames(reference) <- paste0("CellType", 1:n_cell_types)
  
  # True proportions
  true_props <- matrix(
    c(0.5, 0.3, 0.2,
      0.2, 0.5, 0.3,
      0.3, 0.3, 0.4,
      0.6, 0.2, 0.2,
      0.1, 0.6, 0.3),
    nrow = n_samples,
    ncol = n_cell_types,
    byrow = TRUE
  )
  
  # Generate bulk
  bulk <- reference %*% t(true_props) + matrix(rnorm(n_genes * n_samples, sd = 5), nrow = n_genes)
  bulk[bulk < 0] <- 0
  colnames(bulk) <- paste0("Sample", 1:n_samples)
  
  # Run deconvolution
  results <- deconvolve(
    bulk = bulk,
    reference = reference,
    n_cores = 1,
    maxit = 50
  )
  
  # Check output structure
  expect_true(is.list(results))
  expect_true("proportions" %in% names(results))
  expect_true("convergence" %in% names(results))
  
  # Check dimensions
  expect_equal(nrow(results$proportions), n_samples)
  expect_equal(ncol(results$proportions), n_cell_types)
  
  # Check proportions sum to 1
  row_sums <- rowSums(results$proportions)
  expect_true(all(abs(row_sums - 1) < 1e-6))
  
  # Check proportions are non-negative
  expect_true(all(results$proportions >= 0))
})

test_that("deconvolve handles mismatched genes", {
  
  set.seed(42)
  
  # Reference with genes 1-100
  reference <- matrix(rpois(100 * 3, lambda = 100), nrow = 100, ncol = 3)
  rownames(reference) <- paste0("Gene", 1:100)
  colnames(reference) <- paste0("CellType", 1:3)
  
  # Bulk with genes 50-150 (only 50 overlap)
  bulk <- matrix(rpois(100 * 2, lambda = 100), nrow = 100, ncol = 2)
  rownames(bulk) <- paste0("Gene", 50:149)
  colnames(bulk) <- paste0("Sample", 1:2)
  
  # Should work with overlapping genes
  results <- deconvolve(bulk = bulk, reference = reference, n_cores = 1)
  
  expect_equal(nrow(results$proportions), 2)
  expect_equal(ncol(results$proportions), 3)
})

test_that("deconvolve errors with no common genes", {
  
  reference <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(reference) <- c("GeneA", "GeneB", "GeneC")
  colnames(reference) <- paste0("CellType", 1:3)
  
  bulk <- matrix(1:6, nrow = 3, ncol = 2)
  rownames(bulk) <- c("GeneX", "GeneY", "GeneZ")
  colnames(bulk) <- paste0("Sample", 1:2)
  
  expect_error(
    deconvolve(bulk = bulk, reference = reference),
    "No common genes"
  )
})
