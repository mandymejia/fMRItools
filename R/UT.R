#' Upper Triangular Vector to Matrix
#' 
#' Creates a square matrix from a vector which gives the values for the upper 
#'  triangle. By default, the vector is expected to include the diagonal values,
#'  and the values in the upper triangle are copied to the lower triangle so 
#'  that the result is symmetric. 
#'
#' @param x A vector of values for the upper triangular elements of the desired
#'  square matrix.
#' @param diag The values to put on the diagonal, or \code{"x"} (default) if
#'  \code{x} includes the diagonal values.
#' @param LT The values to put on the lower triangle, or \code{"x"} (default)
#'  to use the upper triangular values such that the result is symmetric.
#'  (\code{LT} should not contain the diagonal values, which are provided with
#'  \code{diag} or \code{x}).
#' 
#' 
#' @export
#' @return A square matrix. 
#' 
#' @examples 
#'  UT2mat(seq(10)) # defaults: diag="x", LT="x"
#'  #      [,1] [,2] [,3] [,4]
#'  # [1,]    1    2    4    7
#'  # [2,]    2    3    5    8
#'  # [3,]    4    5    6    9
#'  # [4,]    7    8    9   10
#'  UT2mat(seq(3), LT=8)
#'  #      [,1] [,2]
#'  # [1,]    1    2
#'  # [2,]    8    3
#'  UT2mat(seq(6), diag=0, LT=seq(7,12))
#'  #      [,1] [,2] [,3] [,4]
#'  # [1,]    0    1    2    4
#'  # [2,]    7    0    3    5
#'  # [3,]    8   10    0    6
#'  # [4,]    9   11   12    0
#' UT2mat(rep(-1, 3), diag=c(4,5,6), LT=0)
#'  #      [,1] [,2] [,3]
#'  # [1,]    4   -1   -1
#'  # [2,]    0    5   -1
#'  # [3,]    0    0    6
UT2mat <- function(x, diag="x", LT="x") {
  # Check arguments.
  x_has_diag <- is_1(diag, "character") && (diag=="x")
  if (!x_has_diag) { stopifnot(all(is.numeric(diag) | is.na(diag))) }
  x_for_lt <- is_1(LT, "character") && (LT=="x")
  if (!x_for_lt) { stopifnot(all(is.numeric(LT) | is.na(LT))) }
  if (!is.vector(x) || !is.numeric(x)) stop('`x` must be a numeric vector.')
  
  # Determine V based on M 
  # (solution to quadratic formula since 
  #   M = V*(V+1)/2 (diag=TRUE) or V*(V-1)/2 (diag=FALSE))
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
  if (x_for_lt) {
    mat <- mat + t(mat)
    if (x_has_diag) { diag(mat) <- diag(mat)/2 } # because added twice
  } else {
    mat[lower.tri(mat, diag=FALSE)] <- LT
  }
  if (!x_has_diag) { diag(mat) <- diag }
  
  mat
}

#' Matrix to Upper Triangular Vector
#'
#' @description Returns the vectorized upper triangle of a square matrix
#' @param x A square matrix
#' @param diag Get the diagonal values too? Default: \code{FALSE}
#' 
#' @export
#' @return The vectorized upper triangle of x.  
#'
mat2UT <- function(x, diag=FALSE){
  x[upper.tri(x, diag=diag)]  
}