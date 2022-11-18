#' Transform vector data to an image
#'
#' This fills in parts of a template with values from \code{vec_data}.
#'
#' @param vec_data A V by p matrix, where V is the number of voxels within a
#'   mask and p is the number of vectors to transform into matrix images
#' @param template_image A binary matrix in which V entries are 1 and the rest
#'   of the entries are zero
#'
#' @return A list of masked values from \code{vec_data}
#' 
#' @export
vec2image <- function(vec_data, template_image) {
  each_col <- sapply(split(vec_data, col(vec_data)), function(vd) {
    out <- template_image
    out[out == 1] <- vd
    out[out == 0] <- NA
    return(out)
  }, simplify = F)
  return(each_col)
}

#' Positive skew?
#'
#' Does the vector have a positive skew?
#'
#' @param x The numeric vector for which to calculate the skew. Can also be a matrix,
#'  in which case the skew of each column will be calculated.
#' @return \code{TRUE} if the skew is positive or zero. \code{FALSE} if the skew is negative.
#' @keywords internal
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
#' @param x The ICA results with entries \code{S} and \code{M}
#' @return \code{x} but with positive skew source signals
#' @keywords internal
#'
sign_flip <- function(x){
  stopifnot(is.list(x))
  stopifnot(("S" %in% names(x)) & ("M" %in% names(x)))
  spos <- skew_pos(x$S)
  x$M[,!spos] <- -x$M[,!spos]
  x$S[,!spos] <- -x$S[,!spos]
  x
}

#' Center cols
#'
#' Efficiently center columns of a matrix. (Faster than \code{scale})
#'
#' @param X The data matrix. Its columns will be centered
#' @return The centered data
#' @keywords internal
colCenter <- function(X) {
  X - rep(colMeans(X), rep.int(nrow(X), ncol(X)))
}

#' Unmask matrix data
#'
#' @param dat The data
#' @param mask The mask
#' @param mask_dim Rows, \code{1}, (default) or columns, \code{2}
#' @keywords internal
unmask_mat <- function(dat, mask, mask_dim=1){
  stopifnot(is_posNum(mask_dim))
  stopifnot(mask_dim %in% c(1,2))
  if (mask_dim==1) {
    stopifnot(nrow(dat) == sum(mask))
    mdat <- matrix(NA, nrow=length(mask), ncol=ncol(dat))
    mdat[mask,] <- dat
  } else if (mask_dim==2) {
    stopifnot(ncol(dat) == sum(mask))
    mdat <- matrix(NA, nrow=nrow(dat), ncol=length(mask))
    mdat[,mask] <- dat
  }
  mdat
}

#' Row medians
#' 
#' Use \code{robustbase::rowMedians} package if available, \code{apply} if not.
#'
#' @param x The data matrix
#' @param na.rm Ignore NA values? Default: \code{FALSE}
#' 
#' @return The row medians of \code{x}.
#' 
#' @importFrom stats median
#' @keywords internal
rowMedians2 <- function(x, na.rm = FALSE, ...) {
  if (!requireNamespace("robustbase", quietly=TRUE)) {
    robustbase::rowMedians(x, na.rm=na.rm, ...)
  } else {
    apply(x, 1, median, na.rm=na.rm, ...)
  }
}