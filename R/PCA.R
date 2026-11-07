#' PCA for tall matrix
#'
#' Efficient PCA for a tall matrix (many more rows than columns). Uses the
#'  eigendecomposition of the covariance matrix. The dimensionality of the result can be preset
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
#' @importFrom stats cor
#'
#' @return The decomposition: a list with the eigenvalues \code{d}, the
#'  eigenvectors \code{u}, and (if \code{Vdim > 0}) the unit-norm principal
#'  directions \code{v}.
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
  stopifnot(is_1(Vdim) || identical(Vdim, "Q"))

  nV <- nrow(X) #number of brain locations
  nT <- ncol(X) #number of fMRI volumes (reduce this)
  if(nT > nV) warning('More time points than voxels. Are you sure?')

  if (center) {
    X <- colCenter(X)
  } else {
    TOL <- 1e-8
    if (max(abs(colMeans(X))) > TOL) stop('Columns of X must be centered')
  }

  # Covariance matrix (T x T).
  XtX <- crossprod(X) / (nV - 1)

  # Determine PCA dimensionality.
  if (is.null(Q)) {
    if (!requireNamespace("pesel", quietly = TRUE)) {
      stop("Package \"pesel\" needed to read input data. Please install it", call. = FALSE)
    }
    npc_max <- max(1, min(Q_max, nT - 1))
    Q <- suppressWarnings(pesel::pesel(X, npc.max=npc_max, method='homogenous'))$nPCs
    if (Q == 0) {
      warning("`pesel` estimated zero components; using the number of components with above-average variance explained.")
      comps_var_explained <- eigen(XtX, symmetric=TRUE, only.values=TRUE)$values
      Q <- max(1, sum(comps_var_explained > mean(comps_var_explained)))
    }
  }
  if (Q > nT) { warning("Q > number of columns, so setting Q to ", nT, "."); Q <- nT }
  if (identical(Vdim, "Q")) { Vdim <- Q }
  if (Vdim > Q) { warning("Vdim > Q, so setting Vdim to Q."); Vdim <- Q }

  # Perform dimension reduction.
  out <- tryCatch(
    {
      e <- eigen(XtX, symmetric=TRUE)
      list(d = pmax(e$values, 0), u = e$vectors[, seq_len(Q), drop=FALSE])
    },
    error = function(cond) {
      message(cond)

      # If `eigen` fails, try `svd`.
      z <- tryCatch(svd(XtX, nu=Q, nv=0), error=function(e2) NULL)
      if (!is.null(z)) { return(z) }

      # If `svd` fails:
      # First, try detecting and dropping co-linear columns, and running `svd` again.
      z <- XtX
      zcor <- stats::cor(z)
      zcor[upper.tri(zcor)] <- 0
      diag(zcor) <- 0
      zdrop <- apply(zcor, 2, function(cc) any(abs(cc) > .99))
      if (any(zdrop)) {
        if (sum(!zdrop) < Q) {
          stop(
            sum(zdrop), " co-linear columns detected, but dropping them would ",
            "result in less than Q columns remaining.", call. = FALSE
          )
        }
        cat("Dropping", sum(zdrop), "co-linear columns.\n")
        z <- svd(z[,!zdrop], nu=Q, nv=0)
      } else {
        # If no co-linear columns detected, try `corpcor:fast.svd`.
        if (!requireNamespace("corpcor", quietly = TRUE)) {
          stop(
            "`svd` failed, and the backup routine `corpcor::fast.svd` ",
            "is not available since the Package \"corpcor\" is needed. ",
            "Please install it.", call. = FALSE
          )
        }
        cat("Trying `corpcor::fast.svd`.\n")
        z <- corpcor::fast.svd(z)
        z <- list(d = z$d, u = z$u[, seq_len(Q), drop=FALSE])
      }
      z
    }
  )

  # Compute directions. (out$v would have the directions for XtX, not X.)
  # The directions are scaled to have unit norm.
  vdir <- NULL
  if (Vdim > 0) {
    s <- sqrt(pmax(out$d[seq_len(Vdim)], .Machine$double.eps) * (nV - 1))
    vdir <- (X %*% out$u[, seq_len(Vdim), drop=FALSE]) / rep(s, each=nV)
  }
  out["v"] <- list(vdir)

  out
}