library(DeOPUS)

# ── Helpers ────────────────────────────────────────────────────────────────────

make_data <- function(n_genes = 100, n_cell_types = 3, n_samples = 5, seed = 42) {
  set.seed(seed)
  reference <- matrix(
    rpois(n_genes * n_cell_types, lambda = 100),
    nrow = n_genes, ncol = n_cell_types
  )
  rownames(reference) <- paste0("Gene", seq_len(n_genes))
  colnames(reference) <- paste0("CellType", seq_len(n_cell_types))

  true_props <- matrix(
    runif(n_samples * n_cell_types),
    nrow = n_samples, ncol = n_cell_types
  )
  true_props <- true_props / rowSums(true_props)

  bulk <- reference %*% t(true_props) +
    matrix(rnorm(n_genes * n_samples, sd = 5), nrow = n_genes)
  bulk[bulk < 0] <- 0
  colnames(bulk) <- paste0("Sample", seq_len(n_samples))

  list(bulk = bulk, reference = reference, true_props = true_props)
}

# ── Core functionality ─────────────────────────────────────────────────────────

test_that("deconvolve returns list with proportions and convergence", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_true(is.list(res))
  expect_named(res, c("proportions", "convergence"))
})

test_that("proportions matrix has correct dimensions", {
  d <- make_data(n_genes = 100, n_cell_types = 3, n_samples = 5)
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_equal(nrow(res$proportions), 5)
  expect_equal(ncol(res$proportions), 3)
})

test_that("proportions sum to 1 for every sample", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  row_sums <- rowSums(res$proportions)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("proportions are non-negative", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_true(all(res$proportions >= 0))
})

test_that("output row names match bulk column names", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_equal(rownames(res$proportions), colnames(d$bulk))
})

test_that("output column names match reference column names", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_equal(colnames(res$proportions), colnames(d$reference))
})

test_that("convergence list has one entry per sample", {
  d <- make_data(n_samples = 4)
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_length(res$convergence, 4)
  expect_equal(names(res$convergence), colnames(d$bulk))
})

# ── Parameter variants ─────────────────────────────────────────────────────────

test_that("deconvolve works with power = 1 (MAE loss)", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, power = 1, n_cores = 1, maxit = 50)

  expect_true(all(res$proportions >= 0))
  expect_true(all(abs(rowSums(res$proportions) - 1) < 1e-6))
})

test_that("deconvolve works with maxit = 1 (minimal iterations)", {
  d <- make_data()
  expect_no_error(deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 1))
})

test_that("deconvolve works with custom alpha", {
  d <- make_data()
  res <- deconvolve(d$bulk, d$reference, alpha = 0.1, n_cores = 1, maxit = 50)

  expect_true(all(res$proportions >= 0))
})

test_that("verbose = TRUE emits messages without error", {
  d <- make_data()
  expect_no_error(
    expect_message(
      deconvolve(d$bulk, d$reference, verbose = TRUE, n_cores = 1, maxit = 5)
    )
  )
})

# ── Input coercion ─────────────────────────────────────────────────────────────

test_that("deconvolve accepts data.frame inputs and coerces to matrix", {
  d <- make_data()
  bulk_df <- as.data.frame(d$bulk)
  ref_df  <- as.data.frame(d$reference)

  res <- deconvolve(bulk_df, ref_df, n_cores = 1, maxit = 50)

  expect_true(is.matrix(res$proportions))
  expect_equal(nrow(res$proportions), ncol(d$bulk))
})

# ── Gene overlap handling ──────────────────────────────────────────────────────

test_that("deconvolve handles partial gene overlap", {
  set.seed(1)
  reference <- matrix(rpois(100 * 3, lambda = 100), nrow = 100, ncol = 3)
  rownames(reference) <- paste0("Gene", 1:100)
  colnames(reference) <- paste0("CellType", 1:3)

  bulk <- matrix(rpois(100 * 2, lambda = 100), nrow = 100, ncol = 2)
  rownames(bulk) <- paste0("Gene", 50:149)
  colnames(bulk) <- paste0("Sample", 1:2)

  res <- deconvolve(bulk, reference, n_cores = 1, maxit = 50)

  expect_equal(nrow(res$proportions), 2)
  expect_equal(ncol(res$proportions), 3)
})

test_that("deconvolve errors when no genes overlap", {
  reference <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(reference) <- c("GeneA", "GeneB", "GeneC")
  colnames(reference) <- paste0("CellType", 1:3)

  bulk <- matrix(1:6, nrow = 3, ncol = 2)
  rownames(bulk) <- c("GeneX", "GeneY", "GeneZ")
  colnames(bulk) <- paste0("Sample", 1:2)

  expect_error(deconvolve(bulk, reference), "No common genes")
})

# ── Edge cases ─────────────────────────────────────────────────────────────────

test_that("deconvolve works with a single sample", {
  d <- make_data(n_samples = 1)
  res <- deconvolve(d$bulk, d$reference, n_cores = 1, maxit = 50)

  expect_equal(nrow(res$proportions), 1)
  expect_true(abs(sum(res$proportions) - 1) < 1e-6)
})

test_that("deconvolve works with a single cell type", {
  set.seed(7)
  ref <- matrix(rpois(50, lambda = 100), nrow = 50, ncol = 1)
  rownames(ref) <- paste0("Gene", 1:50)
  colnames(ref) <- "CellType1"

  bulk <- matrix(rpois(50 * 3, lambda = 100), nrow = 50, ncol = 3)
  rownames(bulk) <- paste0("Gene", 1:50)
  colnames(bulk) <- paste0("Sample", 1:3)

  res <- deconvolve(bulk, ref, n_cores = 1, maxit = 50)

  expect_equal(ncol(res$proportions), 1)
  expect_true(all(res$proportions >= 0))
  # With one cell type, proportions are either 1 (converged) or 0 (failed)
  expect_true(all(res$proportions %in% c(0, 1) | abs(res$proportions - 1) < 1e-6))
})

test_that("proportions remain non-negative under heavy noise", {
  set.seed(99)
  n_genes <- 200
  ref <- matrix(rpois(n_genes * 4, lambda = 200), nrow = n_genes, ncol = 4)
  rownames(ref) <- paste0("Gene", seq_len(n_genes))
  colnames(ref) <- paste0("CT", 1:4)

  bulk <- matrix(abs(rnorm(n_genes * 10, mean = 100, sd = 500)), nrow = n_genes)
  rownames(bulk) <- paste0("Gene", seq_len(n_genes))
  colnames(bulk) <- paste0("S", 1:10)

  res <- deconvolve(bulk, ref, n_cores = 1, maxit = 30)

  # Non-negativity must always hold
  expect_true(all(res$proportions >= 0))
  # Row sums are either 1 (converged) or 0 (optimization failed under extreme noise)
  row_sums <- rowSums(res$proportions)
  expect_true(all(abs(row_sums - 1) < 1e-6 | row_sums == 0))
})

# ── Internal utility ───────────────────────────────────────────────────────────

test_that("%||% returns left value when non-NULL", {
  f <- DeOPUS:::`%||%`
  expect_equal(f("a", "b"), "a")
  expect_equal(f(42L, 0L), 42L)
})

test_that("%||% returns right value when left is NULL", {
  f <- DeOPUS:::`%||%`
  expect_equal(f(NULL, "default"), "default")
  expect_null(f(NULL, NULL))
})
