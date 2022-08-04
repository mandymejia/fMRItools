#' Create a mask based on vertices that are invalid
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance. See \code{?create.session} and \code{?is.session}
#'  for more details.
#' @param meanTol,varTol Tolerance for mean and variance of each data location. Locations which
#'  do not meet these thresholds are masked out of the analysis. Defaults: \code{1e-6}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats colVars
#' @return A logical vector indicating valid vertices
#'
#' @export
make_mask <- function(data, meanTol=1e-6, varTol=1e-6, verbose=TRUE){

  # For each BOLD data matrix,
  mask_na <- mask_mean <- mask_var <- rep(TRUE, ncol(data[[1]]$BOLD))
  for (ss in seq(length(data))) {
    dss <- data[[ss]]$BOLD
    # Mark columns with any NA or NaN values for removal.
    dss_na <- is.na(dss) | is.nan(dss)
    mask_na[apply(dss_na, 2, any)] <- FALSE
    # Calculate means and variances of columns, except those with any NA or NaN.
    # Mark columns with mean/var falling under the thresholds for removal.
    mask_mean[mask_na][colMeans(dss[,mask_na,drop=FALSE]) < meanTol] <- FALSE
    mask_var[mask_na][matrixStats::colVars(dss[,mask_na,drop=FALSE]) < varTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    warn_part2 <- if (length(data) > 1) { " in at least one scan.\n" } else { ".\n" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA or NaN values", warn_part2))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean", warn_part2))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance", warn_part2))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}

#' Find nonzero elements in matrix
#'
#' Find nonzero element in a matrix using 2-means clustering
#'
#' @importFrom stats kmeans
#'
#' @param beta_est A vector or matrix of values from which values close to zero
#'  should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values
#'  close to zero are assigned the value of zero. The closeness of a value to
#'  zero is found by performing two-means clustering on the absolute values of
#'  beta_est, and ...
#'
#' @export
#'
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}