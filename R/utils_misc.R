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
    GIFTI2 = "GIFTI2",
    gifti2 = "GIFTI2",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    RDS = "MATRIX",
    data = "MATRIX",
    stop("Unrecognized `format`: \"", format, "\".", call. = FALSE)
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

  if (FORMAT == "GIFTI" || FORMAT == "GIFTI2") {
    if (!requireNamespace("gifti", quietly = TRUE)) {
      stop("Package \"gifti\" needed to work with GIFTI data. Please install it.", call. = FALSE)
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
