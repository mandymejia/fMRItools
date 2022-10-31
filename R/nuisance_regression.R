#' Nuisance regression
#'
#' Performs nuisance regression. The data and design matrix must both be
#'  centered, or an intercept must be included in the design matrix!
#'
#' @param Y The \eqn{T} by \eqn{V} or \eqn{V} by \eqn{T} data.
#' @param design The \eqn{T} by \eqn{Q} matrix of nuisance regressors.
#'
#' @return The data after nuisance regression.
#' 
#' @export
#' 
#' @examples 
#' Y <- matrix(rnorm(700), nrow=100)
#' design <- cbind(seq(100), 1)
#' nuisance_regression(Y, design)
nuisance_regression <- function(Y, design){
  # Z <- design
	# if(nrow(Y) != nrow(Z)) stop('Y and Z must have same number of rows')
 	# invZtZ <- solve(t(Z) %*% Z)  #(Z'Z)^{-1}
	# betahat <- invZtZ %*% t(Z) %*% Y #(Z'Z)^{-1} Z'Y
	# return(Y - Z %*% betahat)

  Y <- as.matrix(Y); design <- as.matrix(design)
  I_m_H <- diag(nrow(design)) - hat_matrix(design)
  if (nrow(Y)==nrow(design)) {
    return(I_m_H %*% Y)
  } else if (ncol(Y)==nrow(design)) {
    return(Y %*% I_m_H)
  } else {
    stop("Y and design are not of compatible dimensions.")
  }
}