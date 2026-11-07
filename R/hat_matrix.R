#' Orthonormal basis for the column space of a design matrix
#'
#' Uses the QR decomposition. Columns of \code{design} that are linear
#'  combinations of other columns are ignored.
#'
#' @param design The \eqn{T} by \eqn{Q} design matrix
#'
#' @return The \eqn{T} by \eqn{r} matrix with orthonormal columns, where
#'  \eqn{r} is the rank of \code{design}
#'
#' @keywords internal
hat_Q <- function(design){
  design <- as.matrix(design)
  if (anyNA(design)) { stop("`design` contains NA values.") }
  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(design)
  r <- qrd$rank
  if (r == 0) { return(matrix(0, nrow(design), 0)) }
  # Deficient columns are pivoted to the end, so the first `r` columns of Q
  #  span the retained columns.
  qr.qy(qrd, diag(1, nrow(design), r))
}

#' Hat matrix
#'
#' Get the hat matrix from a design matrix.
#' 
#' Uses the QR decomposition.
#'
#' @param design The \eqn{T} by \eqn{Q} design matrix
#'
#' @return The \eqn{T} by \eqn{T} hat matrix
#' 
#' @export
#' 
#' @examples 
#' hat_matrix(cbind(seq(100), 1))
hat_matrix <- function(design){
  tcrossprod(hat_Q(design))
}
