#' Dice coefficient
#'
#' @param X,Y Numeric vectors, or matrices with locations along the rows and
#'  networks along the columns. Dice coefficients will be calculated between
#'  each pair of networks between \code{X} and \code{Y}; if \code{Y} is not
#'  provided, \code{X} will be used for \code{Y}.
#' 
#' @return The Dice coefficients in a \eqn{Q_x \times Q_y} matrix, where 
#'  \eqn{Q_x} is the number of columns in \code{X}, \eqn{Q_y} is the number of
#'  columns in \code{Y}, and element \eqn{ij} is the Dice coefficient between 
#'  network \eqn{i} of \code{X} and network \eqn{j} of \code{Y}. 
#'
#' @export
dice_coef <- function(X, Y = NULL) {
  X <- as.matrix(X)

  if (is.null(Y)) {
    Y <- X
  } else {
    Y <- as.matrix(Y)
    stopifnot(nrow(X) == nrow(Y))
  }

  dot_prod <- crossprod(X, Y)
  denom <- outer(colSums(X*X), colSums(Y*Y), "+")
  dice <- ifelse(abs(denom) < 1e-12, 0, (2*dot_prod)/denom)
  dice
}