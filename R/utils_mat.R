#' Transform vector data to image
#'
#' From a \eqn{v \times p} matrix of vectorized data and an \eqn{m \times n} 
#'  image mask with \eqn{v} in-mask locations, create a list of \eqn{p}
#'  \eqn{m \times n} data arrays in which the mask locations are filled
#'  in with the vectorized data values.
#'  
#' Consider using \code{abind::abind} to merge the result into a single
#'  array.
#'
#' @param x \eqn{v \times p} matrix, where \eqn{v} is the number of 
#'  voxels within a mask and \eqn{p} is the number of vectors to transform into
#'  matrix images.
#' @param mask \eqn{m \times n} logical matrix in which \code{v} 
#'  entries are \code{TRUE} and the rest are \code{FALSE}.
#' @param fill_value Out-of-mask value in the output image. Default: 
#'  \code{NA}.
#'
#' @return A list of masked values from \code{x}
#' 
#' @examples
#' x <- unvec_mat(
#'  cbind(seq(3), seq(2,4), seq(3,5)), 
#'  matrix(c(rep(TRUE, 3), FALSE), ncol=2),
#'  0
#' )
#' y <- array(c(1,2,3,0,2,3,4,0,3,4,5,0), dim=c(2,2,3))
#' stopifnot(identical(x[[1]], y[,,1]))
#' stopifnot(identical(x[[2]], y[,,2]))
#' stopifnot(identical(x[[3]], y[,,3]))
#'
#' @export
unvec_mat <- function(x, mask, fill_value=NA) {
  # Check arguments.
  if (length(dim(x))==1) { x <- as.matrix(x) }
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(mask))
  if (is.numeric(mask)) {
    class(mask) <- "logical"
  } else {
    stopifnot(is.logical(mask))
  }
  stopifnot(nrow(x) == sum(mask))
  #stopifnot(length(fill_value)==1)

  sapply(split(x, col(x)), function(vd) {
    out <- mask
    out[mask] <- vd
    out[!mask] <- fill_value
    out
  }, simplify = F)
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

#' Convert data values to percent signal.
#' 
#' Convert data values to percent signal.
#' 
#' @param X a \eqn{T} by \eqn{N} numeric matrix. The columns will be normalized to
#'  percent signal.
#' @param center A function that computes the center of a numeric vector.
#'  Default: \code{median}. Other common options include \code{mean} and 
#'  \code{mode}.
#' @param by Should the center be measured individually for each \code{"column"}
#'  (default), or should the center be the same across \code{"all"} columns?
#' 
#' @return \code{X} with its columns normalized to percent signal. (A value of
#'  85 will represent a -15% signal change.)
#' 
#' @export
pct_sig <- function(X, center=median, by=c("column", "all")){
  stopifnot(is.numeric(X))
  stopifnot(length(dim(X))==2)
  stopifnot(is.function(center))
  by <- match.arg(by, c("column", "all"))

  T_ <- nrow(X); N_ <- ncol(X)
  X <- t(X)

  if (by=="column") {
    m <- apply(X, 1, center)
  } else {
    m <- center(as.numeric(X))
  }

  t(X / m * 100)
}