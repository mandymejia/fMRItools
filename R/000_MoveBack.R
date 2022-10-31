#' Find nonzero elements in matrix
#'
#' Find nonzero element in a matrix using 2-means clustering
#'
#' @importFrom stats kmeans
#'
#' @param beta_est A vector or matrix of values from which values close to zero
#'  should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values
#'  close to zero are assigned the value of zero. The closeness of a value to
#'  zero is found by performing two-means clustering on the absolute values of
#'  beta_est, and ...
#'
#' @export
#'
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}