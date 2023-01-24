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