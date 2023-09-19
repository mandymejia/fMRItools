#' Positive skew?
#'
#' Does the vector have a positive skew?
#'
#' @param x The numeric vector for which to calculate the skew. Can also be a 
#'  matrix, in which case the skew of each column will be calculated.
#' @return \code{TRUE} if the skew is positive or zero. \code{FALSE} if the 
#'  skew is negative.
#' @export
#'
#' @importFrom stats median
skew_pos <- function(x){
  x <- as.matrix(x)
  apply(x, 2, median, na.rm=TRUE) <= colMeans(x, na.rm=TRUE)
}

#' Sign match ICA results
#'
#' Flips all source signal estimates (S) to positive skew
#'
#' @param x The ICA results: a list with entries \code{"S"} and \code{"M"}
#' @return \code{x} but with positive skew source signals
#' @export
#'
sign_flip <- function(x){
  # Check arguments.
  stopifnot(is.list(x))
  stopifnot(("S" %in% names(x)) & ("M" %in% names(x)))
  stopifnot(is.matrix(x$M) && is.matrix(x$S))
  stopifnot(ncol(x$M) == ncol(x$S))

  spos <- skew_pos(x$S)
  x$M[,!spos] <- -x$M[,!spos]
  x$S[,!spos] <- -x$S[,!spos]
  x
}

#' Mode of data vector
#' 
#' Get mode of a data vector. But use the median instead of the mode if all 
#'  data values are unique.
#' 
#' @param x The data vector
#' 
#' @return The mode
#' 
#' @export
#' 
Mode <- function(x) {
  q <- unique(x)
  # Use median instead of the mode if all data values are unique.
  if (length(q) == length(x)) { return(median(x)) }
  q[which.max(tabulate(match(x, q)))]
}

#' Get FORMAT from format
#'
#' @param format the file format
#' @return The file FORMAT
#' @keywords internal
#'
get_FORMAT <- function(format){
  switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    GIFTI = "GIFTI",
    gifti = "GIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    RDS = "MATRIX",
    data = "MATRIX"
  )
}

#' Check required packages for the data format
#'
#' @param FORMAT The data FORMAT
#' @return \code{NULL}, invisibly
#' @keywords internal
check_req_ifti_pkg <- function(FORMAT){
  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
  }

  if (FORMAT == "GIFTI") {
    if (!requireNamespace("gifti", quietly = TRUE)) {
      stop("Package \"gifti\" needed to work with NIFTI data. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
  }

  if (FORMAT == "NIFTI") {
    if (!requireNamespace("RNifti", quietly = TRUE)) {
      stop("Package \"RNifti\" needed to work with NIFTI data. Please install it.", call. = FALSE)
    }
  }

  invisible(NULL)
}

#' Create a mask based on vertices that are invalid
#'
#' @param BOLD A \eqn{V \times T} numeric matrix. Each row is a location.
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{-Inf} for \code{meanTol} (ignore), and \code{1e-6} for
#'  {varTol}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @return A logical vector indicating valid vertices
#'
#' @keywords internal
make_mask <- function(BOLD, meanTol=-Inf, varTol=1e-6, verbose=TRUE){
  stopifnot(is.matrix(BOLD))

  mask_na <- mask_mean <- mask_var <- rep(TRUE, nrow(BOLD))
  # Mark columns with any NA or NaN values for removal.
  mask_na[apply(is.na(BOLD) | is.nan(BOLD), 1, any)] <- FALSE
  # Calculate means and variances of columns, except those with any NA or NaN.
  # Mark columns with mean/var falling under the thresholds for removal.
  mask_mean[mask_na][rowMeans(BOLD[mask_na,,drop=FALSE]) < meanTol] <- FALSE
  if (varTol > 0) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      mask_var[mask_na][matrixStats::rowVars(BOLD[mask_na,,drop=FALSE]) < varTol] <- FALSE
    } else {
      mask_var[mask_na][apply(BOLD[mask_na,,drop=FALSE], 1, var) < varTol] <- FALSE
    }
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA/NaN values.\n"))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean.\n"))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance.\n"))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}