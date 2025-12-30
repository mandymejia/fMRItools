#' Dice overlap
#' 
#' @param X,Y If \code{Y} is \code{NULL}, \code{X} should be a binary matrix, 
#'  e.g. with locations along the rows and networks along the columns. 
#' 
#'  If both are provided, they both must be binary vectors (2 network case). 
#' @return The Dice overlap. 
dice_overlap <- function(X, Y = NULL) {
  if (is.null(Y)) {
      if (!is.matrix(X)) stop("X must be a matrix")
      ux <- unique(as.vector(X))
      if (!all(ux %in% c(0, 1))) stop("Matrix X must be binary")

      XtX   <- crossprod(X)             
      sizes <- diag(XtX)
      denom <- outer(sizes, sizes, "+")
      D <- (2 * XtX) / denom

  } else {
      x <- as.vector(X)
      y <- as.vector(Y)
      if (length(x) != length(y)) stop("x and y must have the same length")
      ux <- unique(c(x, y))
      if (!all(ux %in% c(0, 1))) stop("x and y must be binary")

      inter <- sum(x * y)
      ax <- sum(x)                
      by <- sum(y)                
      D <- (2 * inter) / (ax + by)
  }

  D
}