library(DeOPUS)

# ── sampleData structure ──────────────────────────────────────────────────────

test_that("sampleData loads and is a list", {
  data(sampleData)
  expect_true(is.list(sampleData))
})

test_that("sampleData has required components", {
  data(sampleData)
  expect_true(all(c("bulk", "bulkRatio", "cellTypeExpr", "signature",
                     "singleCellLabels", "singleCellSubjects") %in% names(sampleData)))
})

test_that("sampleData$bulk is a numeric matrix with named dimensions", {
  data(sampleData)
  expect_true(is.matrix(sampleData$bulk))
  expect_true(is.numeric(sampleData$bulk))
  expect_false(is.null(rownames(sampleData$bulk)))
  expect_false(is.null(colnames(sampleData$bulk)))
})

test_that("sampleData$bulk has expected dimensions (11852 genes x 512 samples)", {
  data(sampleData)
  expect_equal(nrow(sampleData$bulk), 11852)
  expect_equal(ncol(sampleData$bulk), 512)
})

test_that("sampleData$cellTypeExpr has expected dimensions (11852 genes x 2 cell types)", {
  data(sampleData)
  expect_equal(nrow(sampleData$cellTypeExpr), 11852)
  expect_equal(ncol(sampleData$cellTypeExpr), 2)
})

test_that("sampleData$bulk and $cellTypeExpr share gene names", {
  data(sampleData)
  common <- intersect(rownames(sampleData$bulk), rownames(sampleData$cellTypeExpr))
  expect_gt(length(common), 0)
})

test_that("sampleData$bulkRatio proportions are between 0 and 1", {
  data(sampleData)
  expect_true(all(sampleData$bulkRatio >= 0))
  expect_true(all(sampleData$bulkRatio <= 1))
})

test_that("sampleData$bulkRatio column sums are approximately 1", {
  data(sampleData)
  col_sums <- colSums(sampleData$bulkRatio)
  expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("sampleData bulk values are non-negative", {
  data(sampleData)
  expect_true(all(sampleData$bulk >= 0))
})

test_that("sampleData bulkRatio rows and cellTypeExpr columns refer to the same cell types", {
  data(sampleData)
  # Same set of names, possibly in different order. Aligning by name is therefore
  # always safe (and necessary, since the order does differ in this dataset).
  expect_setequal(rownames(sampleData$bulkRatio), colnames(sampleData$cellTypeExpr))
})

# ── Integration test: deconvolve on a subset of sampleData ───────────────────
#
# This test verifies the *correlation code logic* — alignment by cell-type name,
# finite outputs, positive sign. It does NOT assert a specific accuracy
# threshold, because the deconvolution algorithm may legitimately change in the
# future. The > 0 assertion exists to catch the misalignment bug where unaligned
# columns produced cor = -1 instead of +1.

test_that("deconvolve runs on sampleData and produces correctly-aligned correlations", {
  data(sampleData)

  set.seed(42)
  idx <- sample(ncol(sampleData$bulk), 10)

  results <- deconvolve(
    bulk      = sampleData$bulk[, idx],
    reference = sampleData$cellTypeExpr,
    n_cores   = 1,
    maxit     = 50
  )

  # Output shape
  expect_equal(nrow(results$proportions), 10)
  expect_equal(ncol(results$proportions), 2)
  expect_true(all(results$proportions >= 0))
  expect_true(all(abs(rowSums(results$proportions) - 1) < 1e-6))

  # Align bulkRatio cell-type ORDER to match results$proportions columns
  ct <- colnames(results$proportions)
  true_props <- t(sampleData$bulkRatio[ct, idx, drop = FALSE])

  cor_values <- sapply(seq_len(nrow(results$proportions)), function(i) {
    cor(results$proportions[i, ], true_props[i, ], method = "pearson")
  })

  # One correlation per sample
  expect_length(cor_values, 10)

  # All finite (catches NA/NaN propagation in the correlation code)
  expect_true(all(is.finite(cor_values)))

  # Sign check: catches the misalignment bug. Naive cor() on reversed cell-type
  # order would give mean = -1 instead of a positive number.
  expect_gt(mean(cor_values), 0)
})
