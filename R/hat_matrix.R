#' Hat matrix
#'
#' Get the hat matrix from a design matrix using the QR decomposition.
#'
#' @param design The \eqn{T} by \eqn{Q} design matrix
#'
#' @return The \eqn{T} by \eqn{T} hat matrix
#' 
#' @export
hat_matrix <- function(design){
  design <- as.matrix(design)
  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(design)
  design <- design[, qrd$pivot[seq_len(qrd$rank)], drop=FALSE]
  qrd <- qr(design)
  Qd <- qr.Q(qrd)
  tcrossprod(Qd)
}