#' Sample bulk RNA-seq deconvolution benchmark dataset
#'
#' A simulated benchmark dataset with known ground-truth cell-type proportions
#' for evaluating and demonstrating the DeOPUS deconvolution method.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{bulk}{Numeric matrix (11852 genes x 512 samples) of simulated bulk
#'     RNA-seq expression values in linear scale. Row names are Ensembl gene IDs
#'     (e.g. ENSG00000000419); column names are Mixture1 through Mixture512.}
#'   \item{bulkRatio}{Numeric matrix (2 cell types x 512 samples) containing
#'     the ground-truth cell-type proportions used to generate \code{bulk}.
#'     Row names are \code{"mature NK T cell"} and \code{"mesenchymal stem cell"}.
#'     Note: row order may differ from the column order of \code{cellTypeExpr},
#'     so align by name before comparing predictions to ground truth.}
#'   \item{cellTypeExpr}{Numeric matrix (11852 genes x 2 cell types) of reference
#'     expression profiles derived from single-cell RNA-seq data.
#'     Column names match (but may be in a different order from) the row names of
#'     \code{bulkRatio}.}
#'   \item{signature}{Numeric matrix (974 signature genes x 2 cell types) of
#'     the most discriminative marker genes for each cell type.}
#'   \item{nCellTypes}{Integer. Number of cell types (2).}
#'   \item{markers}{Named list of marker gene vectors, one entry per cell type.}
#'   \item{sigGenes}{Character vector of signature gene names.}
#'   \item{singleCellExpr}{Sparse matrix (dgCMatrix) of the single-cell RNA-seq
#'     expression data used to derive the reference profiles.}
#'   \item{singleCellLabels}{Character vector of cell-type labels for each cell
#'     in \code{singleCellExpr}.}
#'   \item{singleCellSubjects}{Character vector of donor/subject IDs for each
#'     cell in \code{singleCellExpr}.}
#' }
#'
#' @usage data(sampleData)
#'
#' @examples
#' data(sampleData)
#'
#' # Inspect the dataset
#' dim(sampleData$bulk)          # 11852 x 512
#' dim(sampleData$cellTypeExpr)  # 11852 x 2
#'
#' # Run deconvolution on a small subset for speed
#' set.seed(42)
#' idx <- sample(ncol(sampleData$bulk), 10)
#' results <- deconvolve(
#'   bulk      = sampleData$bulk[, idx],
#'   reference = sampleData$cellTypeExpr,
#'   n_cores   = 1,
#'   maxit     = 50
#' )
#'
#' # Align bulkRatio cell-type ORDER to match results$proportions columns
#' # (bulkRatio rows and cellTypeExpr columns may be in a different order).
#' ct <- colnames(results$proportions)
#' true_props <- t(sampleData$bulkRatio[ct, idx, drop = FALSE])
#'
#' # Per-sample Pearson correlation with ground truth
#' cor_values <- sapply(seq_len(nrow(results$proportions)), function(i) {
#'   cor(results$proportions[i, ], true_props[i, ], method = "pearson")
#' })
#' mean(cor_values, na.rm = TRUE)
#'
#' @source
#' Simulated using negative-binomial expression models with known mixing
#' proportions. See \code{data-raw/prepare_data.R} for the generation script.
"sampleData"
