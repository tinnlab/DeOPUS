#' Deconvolve bulk RNA-seq data using DeOPUS
#'
#' @description
#' Estimates cell-type proportions from bulk RNA-seq data using a reference-based
#' approach with hierarchical shrinkage transformation.
#'
#' @param bulk A numeric matrix of bulk expression data (genes x samples).
#' @param reference A numeric matrix of reference expression profiles (genes x cell types).
#' @param alpha Numeric. Regularization parameter for transformation. Default is 0.01.
#' @param power Numeric. Power for loss function (1 = MAE, 2 = MSE). Default is 2.
#' @param n_cores Integer. Number of cores for parallel processing. Default is 1.
#' @param maxit Integer. Maximum number of optimization iterations. Default is 100.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{proportions}: Matrix of estimated cell-type proportions (samples x cell types)
#'   \item \code{convergence}: List of convergence information for each sample
#' }
#'
#' @details
#' DeOPUS applies a multi-level adaptive transformation combining:
#' \itemize{
#'   \item Local shrinkage: Attenuates influence of high-variance genes
#'   \item Global shrinkage: Provides overall regularization
#'   \item Power transformation: Stabilizes variance
#'   \item Quantile normalization: Ensures robust profile comparison
#' }
#'
#' @examples
#' # Load the built-in simulated benchmark dataset
#' data(simulated)
#'
#' # Run deconvolution on a small subset for speed
#' set.seed(42)
#' idx <- sample(ncol(simulated$bulk), 10)
#' results <- deconvolve(
#'   bulk      = simulated$bulk[, idx],
#'   reference = simulated$cellTypeExpr,
#'   n_cores   = 1,
#'   maxit     = 50
#' )
#'
#' # View estimated proportions (samples x cell types)
#' head(results$proportions)
#'
#' @importFrom parallel mclapply
#' @importFrom stats optim qnorm
#' @export
deconvolve <- function(bulk,
                       reference,
                       alpha = 0.01,
                       power = 2,
                       n_cores = 1,
                       maxit = 100,
                       verbose = FALSE) {

  # Input validation
  if (!is.matrix(bulk)) {
    bulk <- as.matrix(bulk)
  }
  if (!is.matrix(reference)) {
    reference <- as.matrix(reference)
  }

  # Check for common genes
  common_genes <- intersect(rownames(bulk), rownames(reference))
  if (length(common_genes) == 0) {
    stop("No common genes found between bulk and reference matrices.")
  }

  if (verbose) {
    message(sprintf("Using %d common genes for deconvolution", length(common_genes)))
  }

  # Subset to common genes
  bulk <- bulk[common_genes, , drop = FALSE]
  reference <- reference[common_genes, , drop = FALSE]

  # Get dimensions
  n_samples <- ncol(bulk)
  n_cell_types <- ncol(reference)
  cell_type_names <- colnames(reference)
  sample_names <- colnames(bulk)

  # Define transformation function
  transform_func <- function(q, alpha = 0.01) {
    q[q < 0] <- 0

    # Normalize by total counts
    scale_factor <- mean(q) * 100
    if (scale_factor == 0) scale_factor <- 1

    # Shrinkage parameters
    tau <- 0.05     # Local shrinkage parameter
    phi <- 0.8      # Global shrinkage intensity
    omega <- 1.5    # Shape modulation factor
    delta <- 0.3    # Adaptive regularization strength

    # Dispersion and shape parameters
    kappa <- 0.1
    mu <- 2.0

    # Hierarchical shrinkage adjustment
    local_shrink <- tau / (tau + q / scale_factor)
    global_shrink <- (1 - phi) * exp(-delta * q / scale_factor) + phi

    # Modified power transformation with adaptive hierarchical shrinkage
    z <- ((q / scale_factor + kappa)^(mu * omega * global_shrink) - kappa^mu) / (mu * global_shrink)

    # Apply spline-based normalization
    ranks <- rank(z) / length(z)
    quantile_norm <- stats::qnorm(pmax(0.001, pmin(0.999, ranks)))

    # Enhanced blending transformation with multi-level shrinkage
    signal_weight <- pmin(1, q / (10 + alpha * q)) * (1 - local_shrink)

    # Incorporate prior information shrinkage
    transformed <- signal_weight * z + (1 - signal_weight) *
      global_shrink *
      quantile_norm

    # Apply edge-preserving smoothing via shrinkage
    transformed <- transformed / (1 + delta * abs(transformed))

    # Scale to approximate log2 scale with adaptive shrinkage
    transformed <- transformed / max(abs(transformed), na.rm = TRUE) * log2(max(q) + 1)

    return(transformed)
  }

  # Define loss function
  loss_function <- function(p, b, reference, alpha, power, transform_func) {
    b_transformed <- transform_func(b, alpha)
    fitted <- as.vector(reference %*% p)
    fitted_transformed <- transform_func(fitted, alpha)
    return(sum(abs(b_transformed - fitted_transformed)^power))
  }

  # Run deconvolution for each sample
  if (verbose) {
    message(sprintf("Deconvolving %d samples using %d cores...", n_samples, n_cores))
  }

  results <- parallel::mclapply(seq_len(n_samples), mc.cores = n_cores, mc.preschedule = FALSE, function(i) {
    tryCatch({
      optim_result <- stats::optim(
        par = rep(1 / n_cell_types, n_cell_types),
        fn = loss_function,
        b = bulk[, i],
        reference = reference,
        alpha = alpha,
        power = power,
        transform_func = transform_func,
        method = "L-BFGS-B",
        lower = 0,
        upper = 100,
        control = list(maxit = maxit)
      )
      list(
        proportions = optim_result$par,
        convergence = optim_result$convergence,
        value = optim_result$value
      )
    }, error = function(e) {
      list(
        proportions = rep(0, n_cell_types),
        convergence = -1,
        value = NA,
        error = conditionMessage(e)
      )
    })
  })

  # Extract proportions matrix
  proportions <- do.call(rbind, lapply(results, function(x) x$proportions))

  # Normalize to sum to 1
  row_sums <- rowSums(proportions)
  row_sums[row_sums == 0] <- 1
  proportions <- proportions / row_sums

  # Handle NAs
  proportions[is.na(proportions)] <- 0

  # Set names
  colnames(proportions) <- cell_type_names
  rownames(proportions) <- sample_names

  # Extract convergence info
  convergence <- lapply(results, function(x) {
    list(
      convergence = x$convergence,
      value = x$value,
      error = x$error %||% NULL
    )
  })
  names(convergence) <- sample_names

  if (verbose) {
    n_converged <- sum(sapply(results, function(x) x$convergence == 0))
    message(sprintf("Deconvolution complete. %d/%d samples converged.", n_converged, n_samples))
  }

  return(list(
    proportions = proportions,
    convergence = convergence
  ))
}

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
