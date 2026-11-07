#' Create a mask based on vertices that are invalid
#'
#' @param BOLD A \eqn{V \times T} numeric matrix. Each row is a location.
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{-Inf} for \code{meanTol} (ignore), and \code{1e-6} for
#'  \code{varTol}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats rowVars rowAnyNAs
#' @return A logical vector indicating valid vertices
#'
#' @export
mask_BOLD <- function(BOLD, meanTol=-Inf, varTol=1e-6, verbose=TRUE){
  stopifnot(is.matrix(BOLD))

  # Mark locations with any NA or NaN values for removal.
  mask_na <- !matrixStats::rowAnyNAs(BOLD)
  mask_mean <- mask_var <- rep(TRUE, nrow(BOLD))
  # Calculate means and variances of locations, except those with any NA or NaN.
  # Mark locations with mean/var falling under the thresholds for removal.
  if (any(mask_na)) {
    BOLD_valid <- if (all(mask_na)) { BOLD } else { BOLD[mask_na,,drop=FALSE] }
    if (is.finite(meanTol)) {
      mask_mean[mask_na] <- rowMeans(BOLD_valid) >= meanTol
    }
    if (varTol > 0) {
      v <- matrixStats::rowVars(BOLD_valid)
      mask_var[mask_na] <- !is.na(v) & v >= varTol
    }
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    n_na <- sum(!mask_na)
    # Do not include NA locations in count.
    n_mean <- sum(!mask_mean)
    # Do not include NA or low-mean locations in count.
    n_var <- sum(!mask_var & mask_mean)
    warn_part1 <- if (n_na > 0) { "additional locations" } else { "locations" }
    if (n_na > 0) {
      cat(n_na, paste0("locations removed due to NA/NaN values.  "))
    }
    if (n_mean > 0) {
      cat(n_mean, warn_part1, paste0("removed due to low mean.  "))
    }
    if (n_var > 0) {
      cat(n_var, warn_part1, paste0("removed due to low variance.  "))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}
