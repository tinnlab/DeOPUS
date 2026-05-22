## Script that generated data/sampleData.rda
##
## Run this script from the package root to regenerate the package dataset
## from the original source file.

# ── sampleData ────────────────────────────────────────────────────────────────
#
# The source file was a pre-generated simulation containing:
#   11852 genes x 512 simulated bulk RNA-seq samples
#   2 cell types (mature NK T cell, mesenchymal stem cell)
#   Matching ground-truth proportions (bulkRatio)

# Load original source (adjust path as needed)
sampleData <- readRDS("path/to/source_simulation.rds")

# Save as compressed .rda for use as package data
usethis::use_data(sampleData, overwrite = TRUE, compress = "xz")
