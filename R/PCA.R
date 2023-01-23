#' PCA for tall matrix
#' 
#' Efficient PCA for a tall matrix (many more rows than columns). Uses the SVD
#'  of the covariance matrix. The dimensionality of the result can be preset
#'  with \code{Q} or estimated with PESEL. 
#' 
#' @param X The tall numeric matrix for which to compute the PCA. For fMRI data, 
#'  \code{X} should be \code{V} brain locations by \code{T} timepoints.
#' @param center Center the columns of \code{X}? Default: \code{TRUE}. Set to
#'  \code{FALSE} if already centered. Centered data is required to compute PCA. 
#' @param Q Number of latent dimensions to estimate. If \code{NULL} (default), 
#'  estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic 
#'  dimensionality selection with PESEL. Default: \code{100}.
#' @param Vdim Number of principal directions to obtain. Default: \code{0}. Can
#'  also be \code{"Q"} to set equal to the value of \code{Q}. Note that setting
#'  this value less than \code{Q} does not speed up computation time, but does
#'  save on memory. Note that the directions will be with respect to \code{X},
#'  not its covariance matrix.
#' 
#' @export
#' 
#' @return The SVD decomposition
#' 
#' @examples 
#' U <- matrix(rnorm(900), nrow=300, ncol=3)
#' V <- matrix(rnorm(15), nrow=3, ncol=5)
#' PCA(U %*% V)
PCA <- function(X, center=TRUE, Q=NULL, Q_max=100, Vdim=0) {
  # Check arguments.
  stopifnot(is.matrix(X) && is.numeric(X))
  stopifnot(is_1(center, "logical"))
  stopifnot(is.null(Q) || is_1(Q))
  stopifnot(is_1(Q_max))
  stopifnot(is_1(Vdim))

  nV <- nrow(X) #number of brain locations
  nT <- ncol(X) #number of fMRI volumes (reduce this)
  if(nT > nV) warning('More time points than voxels. Are you sure?')

  if (center) {
    X <- colCenter(X)
  } else {
    TOL <- 1e-8
    if (max(abs(colMeans(X))) > TOL) stop('Columns of X must be centered')
  }
  
  # Determine PCA dimensionality.
  if(is.null(Q)){
    if (!requireNamespace("pesel", quietly = TRUE)) {
      stop("Package \"pesel\" needed to read input data. Please install it", call. = FALSE)
    }
    Q <- suppressWarnings(pesel::pesel(X, npc.max=Q_max, method='homogenous'))$nPCs
  }
  if (Vdim == "Q") { Vdim <- Q }
  if (Vdim > Q) { warning("Vdim > Q, so setting Vdim to Q."); Vdim <- Q }

  # Perform dimension reduction.
  out <- svd(crossprod(X) / (nV-1), nu=Q, nv=0)

  # Compute directions. (out$v would have the directions for XtX, not X.)
  if (Vdim > 0) {
    out$v <- X %*% out$u[,seq(Vdim)] %*% diag(1/out$d[seq(Vdim)])
  } else {
    out["v"] <- list(NULL)
  }

  out
}