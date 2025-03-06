#' Upper Triangular Vector to Matrix
#' @description  Returns the symmatric square matrix from a vector containing 
#' the upper triangular elements
#'
#' @param x A vector containing the upper triangular elements of a square, 
#' symmetric matrix.
#' @param diag A scalar value to use for the diagonal values of the matrix, or
#'  \code{"x"} if \code{x} includes the diagonal values. Default: \code{1}.
#' @param LT Change from \code{TRUE} (default) to \code{FALSE} to set lower
#'  triangle values to zero.
#' 
#' @export
#' @return If \code{LT}, a symmetric matrix with the values of \code{x} in the 
#'  upper and lower triangles and the value \code{diag} on the diagonal. If
#'  \code{!LT}, the lower triangle values will be zero instead.
#' 
UT2mat <- function(x, diag=1, LT=TRUE) {
  stopifnot(length(diag)==1)
  x_has_diag <- is.character(diag) && (diag=="x")
  if (!x_has_diag) { stopifnot(is.numeric(diag) | is.na(diag)) }

  if(!is.vector(x) | !is.numeric(x)) stop('`x` must be a numeric vector.')
  
  #determine V based on M (solution to quadratic formula since M = V*(V+1)/2 (diag=TRUE) or V*(V-1)/2 (diag=FALSE))
  M <- length(x)
  V <- if (x_has_diag) { (-1+sqrt(8*M+1))/2 } else { (1+sqrt(8*M+1))/2 }
  if (round(V) != V) {
    if (x_has_diag) {
      stop('Length of x not equal to V(V+1)/2, for some integer V.')
    } else {
      stop('Length of x not equal to V(V-1)/2, for some integer V.')
    }
  }
  mat <- matrix(0, nrow=V, ncol=V)
  mat[upper.tri(mat, diag=x_has_diag)] <- x
  if (LT) {
    mat <- mat + t(mat)
    if (x_has_diag) { diag(mat) <- diag(mat)/2 }
  }
  if (!x_has_diag) { diag(mat) <- diag }
  
  mat
}

#' Matrix to Upper Triangular Vector
#'
#' @description Returns the vectorized upper triangle of a square matrix
#' @param x A square matrix
#' 
#' @export
#' @return The vectorized upper triangle of x.  
#'
mat2UT <- function(x){
  x[upper.tri(x)]  
}