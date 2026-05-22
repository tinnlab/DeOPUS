library(DeOPUS)

# ── simulated dataset ──────────────────────────────────────────────────────────

test_that("simulated dataset loads and is a list", {
  data(simulated)
  expect_true(is.list(simulated))
})

test_that("simulated dataset has required components", {
  data(simulated)
  expect_true(all(c("bulk", "bulkRatio", "cellTypeExpr", "signature",
                     "singleCellLabels", "singleCellSubjects") %in% names(simulated)))
})

test_that("simulated$bulk is a numeric matrix with named dimensions", {
  data(simulated)
  expect_true(is.matrix(simulated$bulk))
  expect_true(is.numeric(simulated$bulk))
  expect_false(is.null(rownames(simulated$bulk)))
  expect_false(is.null(colnames(simulated$bulk)))
})

test_that("simulated$bulk has expected dimensions (11852 genes x 512 samples)", {
  data(simulated)
  expect_equal(nrow(simulated$bulk), 11852)
  expect_equal(ncol(simulated$bulk), 512)
})

test_that("simulated$cellTypeExpr has expected dimensions (11852 genes x 2 cell types)", {
  data(simulated)
  expect_equal(nrow(simulated$cellTypeExpr), 11852)
  expect_equal(ncol(simulated$cellTypeExpr), 2)
})

test_that("simulated$bulk and $cellTypeExpr share gene names", {
  data(simulated)
  common <- intersect(rownames(simulated$bulk), rownames(simulated$cellTypeExpr))
  expect_gt(length(common), 0)
})

test_that("simulated$bulkRatio proportions are between 0 and 1", {
  data(simulated)
  expect_true(all(simulated$bulkRatio >= 0))
  expect_true(all(simulated$bulkRatio <= 1))
})

test_that("simulated$bulkRatio column sums are approximately 1", {
  data(simulated)
  col_sums <- colSums(simulated$bulkRatio)
  expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("simulated bulk values are non-negative", {
  data(simulated)
  expect_true(all(simulated$bulk >= 0))
})

test_that("deconvolve runs successfully on a subset of simulated data", {
  data(simulated)
  set.seed(42)
  idx <- sample(ncol(simulated$bulk), 5)

  res <- deconvolve(
    bulk      = simulated$bulk[, idx],
    reference = simulated$cellTypeExpr,
    n_cores   = 1,
    maxit     = 30
  )

  expect_true(is.list(res))
  expect_equal(nrow(res$proportions), 5)
  expect_equal(ncol(res$proportions), 2)
  expect_true(all(res$proportions >= 0))
  expect_true(all(abs(rowSums(res$proportions) - 1) < 1e-6))
})

test_that("deconvolve results on simulated data have correct cell-type names", {
  data(simulated)
  set.seed(1)
  idx <- sample(ncol(simulated$bulk), 3)

  res <- deconvolve(
    bulk      = simulated$bulk[, idx],
    reference = simulated$cellTypeExpr,
    n_cores   = 1,
    maxit     = 30
  )

  expect_equal(sort(colnames(res$proportions)),
               sort(colnames(simulated$cellTypeExpr)))
})

# ── GSE77343 dataset (inst/extdata) ───────────────────────────────────────────

test_that("GSE77343.rds is accessible via system.file", {
  path <- system.file("extdata", "GSE77343.rds", package = "DeOPUS")
  expect_true(nchar(path) > 0)
  expect_true(file.exists(path))
})

test_that("GSE77343 loads as a list with required components", {
  path <- system.file("extdata", "GSE77343.rds", package = "DeOPUS")
  gse  <- readRDS(path)

  expect_true(is.list(gse))
  expect_true(all(c("bulk", "cellTypeExpr", "bulkRatio") %in% names(gse)))
})

test_that("GSE77343$bulk is a numeric matrix with named dimensions", {
  path <- system.file("extdata", "GSE77343.rds", package = "DeOPUS")
  gse  <- readRDS(path)

  expect_true(is.matrix(gse$bulk))
  expect_true(is.numeric(gse$bulk))
  expect_false(is.null(rownames(gse$bulk)))
  expect_false(is.null(colnames(gse$bulk)))
})

test_that("GSE77343$bulk and $cellTypeExpr share gene names", {
  path <- system.file("extdata", "GSE77343.rds", package = "DeOPUS")
  gse  <- readRDS(path)

  common <- intersect(rownames(gse$bulk), rownames(gse$cellTypeExpr))
  expect_gt(length(common), 0)
})
