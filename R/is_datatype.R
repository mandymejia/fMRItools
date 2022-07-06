#' Is this an integer?
#' 
#' @param x The putative integer
#' @param nneg Require \code{x>=0} (non-negative) too?
#' @return Logical indicating whether \code{x} is an integer
#' 
#' @keywords internal
is_integer <- function(x, nneg=FALSE){
  y <- FALSE
  if (length(x)==1 && is.numeric(x)) {
    if (x%%1==0) {
      if (x>=0 || !nneg) { y <- TRUE }
    }
  } 
  y
}

#' Is this numeric vector constant?
#' 
#' @param x The numeric vector
#' @param TOL minimum range of \code{x} to be considered non-constant.
#'  Default: \code{1e-8}
#' 
#' @return Is \code{x} constant? 
#' 
#' @keywords internal
is_constant <- function(x, TOL=1e-8) {
  abs(max(x) - min(x)) < TOL
}