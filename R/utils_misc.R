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

#' Center matrix columns
#'
#' Efficiently center columns of a matrix. (Faster than \code{base::scale}.)
#'
#' @param X The data matrix. Its columns will be centered.
#' @return The centered data
#' @export
colCenter <- function(X) {
  X - rep(colMeans(X), rep.int(nrow(X), ncol(X)))
}

#' Unmask matrix data
#' 
#' Insert empty rows or columns to a matrix. For example, medial wall vertices 
#'  can be added back to the cortex data matrix. 
#'
#' @param x The data matrix to unmask.
#' @param mask The logical mask: the number of \code{TRUE} values should match
#'  the size of the (\code{mask_dim})th dimension in \code{dat}.
#' @param mask_dim Rows, \code{1} (default), or columns, \code{2}.
#' @param fill The fill value for the inserted rows/columns. Default: \code{NA}.
#' 
#' @return The unmasked matrix.
#' @export
unmask_mat <- function(x, mask, mask_dim=1, fill=NA){
  # Argument checks.
  stopifnot(is.matrix(x))
  stopifnot(is.logical(mask))
  stopifnot(is_posNum(mask_dim))
  stopifnot(mask_dim %in% c(1,2))
  stopifnot(length(fill) == 1)

  # Unmask.
  if (mask_dim==1) {
    stopifnot(nrow(x) == sum(mask))
    mdat <- matrix(fill, nrow=length(mask), ncol=ncol(x))
    mdat[mask,] <- x
  } else if (mask_dim==2) {
    stopifnot(ncol(x) == sum(mask))
    mdat <- matrix(fill, nrow=nrow(x), ncol=length(mask))
    mdat[,mask] <- x
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