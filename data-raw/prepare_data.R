## Script that generated data/simulated.rda and inst/extdata/GSE77343.rds
##
## Run this script from the package root to regenerate the package datasets
## from the original source files.

# ── Simulated dataset ────────────────────────────────────────────────────────
#
# The original source file was provided as a pre-generated simulation:
#   7c3e2a286fda1722fe591fc527ded32c.rds (17 MB)
#
# It contains 11852 genes x 512 simulated bulk RNA-seq samples with 2 cell
# types (mature NK T cell, mesenchymal stem cell) and matching ground-truth
# proportions (bulkRatio).

# Load original source (adjust path as needed)
simulated <- readRDS("path/to/7c3e2a286fda1722fe591fc527ded32c.rds")

# Save as compressed .rda for use as package data
usethis::use_data(simulated, overwrite = TRUE, compress = "xz")

# ── GSE77343 real dataset ─────────────────────────────────────────────────────
#
# Real bulk RNA-seq data from GEO accession GSE77343 (Linsley et al. 2014).
# 17275 genes x 197 samples with 15 reference cell types.
# Original file: GSE77343.rds (88 MB)
#
# Due to its size, this dataset is stored in inst/extdata/ rather than data/,
# and must be loaded explicitly:
#   readRDS(system.file("extdata", "GSE77343.rds", package = "DeOPUS"))

file.copy(
  from = "path/to/GSE77343.rds",
  to   = file.path("inst", "extdata", "GSE77343.rds"),
  overwrite = TRUE
)
