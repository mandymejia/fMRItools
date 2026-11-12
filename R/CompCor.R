#' Construct a CompCor noise object
#'
#' Common return type for the noise PCs, shared by \code{\link{CompCor}} (its
#'  \code{"noise"} entry) and \code{\link{CompCor_HCP_calc}}.
#'
#' @param PCs,var,ROI_noise The noise PC scores, their variances, and the noise
#'  ROI masks (after erosion), each a list with one entry per noise ROI.
#'
#' @return A list with entries \code{"PCs"}, \code{"var"} and \code{"ROI_noise"},
#'  of class \code{"CompCor_noise"}.
#'
#' @keywords internal
new_CompCor_noise <- function(PCs, var, ROI_noise) {
  structure(
    list(PCs=PCs, var=var, ROI_noise=ROI_noise),
    class="CompCor_noise"
  )
}

#' CompCor: get noise PCs
#'
#' Computes the noise principal components (PCs) for aCompCor from the data of
#'  each noise ROI. This is the "compute" step shared by \code{\link{CompCor}}
#'  and \code{\link{CompCor_HCP_calc}}. The "apply" step is
#'  \code{\link{CompCor_clean}}.
#'
#' If \code{prelim_nuisance} is provided, it is regressed from every noise ROI
#'  first (the orthonormal basis of its column space is computed only once, and
#'  rank-deficient designs are handled). Then, for each ROI, the voxels are
#'  centered and scaled (see \code{center} and \code{scale}) and the PCs are
#'  obtained with the singular value decomposition.
#'
#' Note that centering is by the median, which is not the mean. If
#'  \code{prelim_nuisance} includes an intercept, the preliminary regression
#'  makes each voxel's mean zero, but subtracting the median afterwards can add
#'  a small constant back. The PCs are therefore orthogonal to the non-constant
#'  part of \code{prelim_nuisance}, but may retain a small component along the
#'  constant. (With \code{center=FALSE} they are orthogonal to all of
#'  \code{prelim_nuisance}.)
#'
#' @param X_noise A list of the noise ROIs data. Each entry is a \eqn{T} by
#'  \eqn{V} matrix (timepoints by voxels), with \eqn{V} possibly differing
#'  across ROIs.
#' @param center,scale Center the columns of the noise ROI data by their
#'  medians, and scale by their MADs? Default: \code{TRUE} for both.
#' @param noise_nPC The number of PCs to compute for each noise ROI: a single
#'  number, or a list or vector with one entry per ROI. Values between 0 and 1
#'  represent the minimum proportion of variance explained by the PCs used for
#'  each noise ROI. Default: \code{5}.
#' @param prelim_nuisance Nuisance signals to regress from the noise ROIs before
#'  computing the PCs. A \eqn{T} by \eqn{P} numeric matrix. Default:
#'  \code{NULL} (do not do any preliminary nuisance regression).
#'
#' @return A list with entries \code{"PCs"} (the noise PC scores), \code{"var"}
#'  (the variance of each PC) and \code{"vartotal"} (the total variance), each a
#'  list with one entry per noise ROI.
#'
#' @keywords internal
noise_PCs <- function(
  X_noise, center=TRUE, scale=TRUE, noise_nPC=5, prelim_nuisance=NULL){
  TOL <- 1e-8

  N <- length(X_noise)
  PCs <- vector("list", N); names(PCs) <- names(X_noise)
  noise_var <- vector("list", N); names(noise_var) <- names(X_noise)
  noise_vartotal <- vector("list", N); names(noise_vartotal) <- names(X_noise)

  if (length(noise_nPC) == 1) {
    noise_nPC <- rep(noise_nPC, N)
  }
  noise_nPC <- as.list(noise_nPC)
  if (is.null(names(noise_nPC))) {
    if (length(noise_nPC) < N) {
      stop("The length of `noise_nPC` is less than the number of noise ROIs.")
    }
    noise_nPC <- noise_nPC[seq_len(N)]
    names(noise_nPC) <- names(X_noise)
  } else {
    if (!all(names(X_noise) %in% names(noise_nPC))) {
      stop("The names of `noise_nPC` do not match those of the noise ROIs.")
    }
    noise_nPC <- noise_nPC[names(X_noise)]
  }

  # Preliminary nuisance regression: project out the column space of the
  #   design. `hat_Q` gives an orthonormal basis for it, which is cheaper than
  #   the T x T hat matrix and is only computed once for all ROIs.
  Qd <- NULL
  if (!is.null(prelim_nuisance)) {
    prelim_nuisance <- as.matrix(prelim_nuisance)
    Qd <- hat_Q(prelim_nuisance)
  }

  for (ii in seq_len(N)) {
    if (ncol(X_noise[[ii]])==0) { next }

    if (!is.null(Qd)) {
      if (nrow(X_noise[[ii]]) != nrow(prelim_nuisance)) {
        stop(
          "`prelim_nuisance` has ", nrow(prelim_nuisance), " rows, but the ",
          "noise ROI data has ", nrow(X_noise[[ii]]), " timepoints."
        )
      }
      X_noise[[ii]] <- X_noise[[ii]] - Qd %*% crossprod(Qd, X_noise[[ii]])
    }

    # Transpose.
    X_noise[[ii]] <- t(X_noise[[ii]])
    #	Center.
    if (center) { X_noise[[ii]] <- X_noise[[ii]] - c(rowMedians2(X_noise[[ii]], na.rm=TRUE)) }
    # Compute MADs.
    mad <- 1.4826 * rowMedians2(abs(X_noise[[ii]]), na.rm=TRUE)
    X_constant <- mad < TOL
    if (any(X_constant)) {
      if (all(X_constant)) {
      stop(paste0("All data locations in noise ROI ", ii, " are zero-variance.\n"))
      } else {
        warning(paste0("Warning: removing ", sum(X_constant),
        " constant data locations (out of ", length(X_constant),
        ") in noise ROI ", ii,
        ".\n"))
      }
    }
    mad <- mad[!X_constant]; X_noise[[ii]] <- X_noise[[ii]][!X_constant,]
    # Scale.
    if (scale) { X_noise[[ii]] <- X_noise[[ii]]/c(mad) }
    # Revert transpose.
    X_noise[[ii]] <- t(X_noise[[ii]])
    if (ncol(X_noise[[ii]])==0) { next }

    # Compute the PC scores.
    x <- svd(tcrossprod(X_noise[[ii]]))
    noise_var[[ii]] <- x$d
    noise_vartotal[[ii]] <- sum(noise_var[[ii]])
    if (noise_nPC[[ii]] < 1) {
      # Use enough PCs to explain the desired proportion of variance.
      noise_nPC[[ii]] <- min(
        length(x$d),
        sum(cumsum(noise_var[[ii]]) < noise_vartotal[[ii]]*noise_nPC[[ii]]) + 1
      )
    }
    PCs[[ii]] <- x$u[,seq(noise_nPC[[ii]]),drop=FALSE]
    noise_var[[ii]] <- noise_var[[ii]][seq(noise_nPC[[ii]])]
  }

  list(PCs=PCs, var=noise_var, vartotal=noise_vartotal)
}

#' CompCor: clean data with noise PCs
#'
#' Removes the noise PCs from data: the "apply" step of aCompCor. The data is
#'  regressed on a design matrix with the noise PCs, any additional regressors
#'  (DCT bases and/or \code{nuisance}), and an intercept term. The PCs can be
#'  obtained with \code{\link{CompCor}} or \code{\link{CompCor_HCP_calc}}.
#'  \code{\link{CompCor}} calls this function, and so does
#'  \code{\link{CompCor_HCP_apply}}, which also handles CIFTI data.
#'
#' @param X The \eqn{T \times V} or \eqn{V \times T} data matrix, where \eqn{T}
#'  is the number of timepoints. (If square, it is assumed to be \eqn{T \times
#'  V}.) The data does not need to be centered, since an intercept is always
#'  included in the regression.
#' @param noise The noise PCs: a list with an entry named \code{"PCs"} which is
#'  a list of \eqn{T} by \eqn{Q_i} matrices, one for each noise ROI. For
#'  example, the result of \code{\link{CompCor_HCP_calc}} or the \code{"noise"}
#'  entry of the result of \code{\link{CompCor}}.
#' @param DCT Add DCT bases to the nuisance regression? Use an integer to
#'  indicate the number of cosine bases. Use \code{0} (default) to forgo
#'  detrending.
#' @param nuisance Nuisance signals to regress from each data column in addition
#'  to the noise ROI PCs. Should be a \eqn{T} by \eqn{N} numeric matrix where
#'  \eqn{N} represents the number of nuisance signals. To not perform any
#'  additional nuisance regression set this argument to \code{NULL}, \code{0}, or
#'  \code{FALSE}. Default: \code{NULL}.
#'
#' @return The cleaned data, in the same orientation as \code{X}.
#'
#' @keywords internal
#' @seealso CompCor, CompCor_HCP_calc, CompCor_HCP_apply
CompCor_clean <- function(X, noise, DCT=0, nuisance=NULL){
  X <- as.matrix(X)
  stopifnot(is.list(noise), !is.null(noise$PCs))
  stopifnot(is.numeric(DCT), length(DCT)==1, DCT >= 0)

  # Make the design matrix.
  design <- do.call(cbind, noise$PCs)
  if (is.null(design)) { stop("`noise` does not contain any PCs.") }
  T_ <- nrow(design)
  if (!(nrow(X) == T_ || ncol(X) == T_)) {
    stop(
      "The noise PCs have ", T_, " timepoints, but `X` has dimensions ",
      nrow(X), " by ", ncol(X), "."
    )
  }
  if (DCT > 0) { design <- cbind(design, dct_bases(T_, DCT)) }
  if (!(is.null(nuisance) || isFALSE(nuisance) || identical(nuisance, 0))) {
    if (nrow(nuisance) != T_) {
      stop("`nuisance` does not have the same number of timepoints as the noise PCs.")
    }
    nuisance <- validate_design_mat(nuisance, T_)
    design <- cbind(design, nuisance)
  }
  design <- validate_design_mat(cbind(1, design), T_)

  nuisance_regression(X, design)
}

#' Anatomical CompCor
#'
#' The aCompCor algorithm for denoising fMRI data using noise ROIs data
#'
#' This function performs two steps, each of which is available as a separate
#'  function. First, the principal components (PCs) of each noise region of
#'  interest (ROI) are calculated. For each ROI, voxels are centered and scaled
#'  (can be disabled with the arguments \code{center} and \code{scale}), and then
#'  the PCs are calculated via the singular value decomposition.
#'
#' Next, aCompCor is performed to remove the shared variation between the noise ROI
#'  PCs and each location in the data. This is accomplished by a nuisance regression
#'  using a design matrix with the noise ROI PCs, any additional regressors specified
#'  by \code{nuisance}, and an intercept term, which is done by
#'  \code{\link{CompCor_clean}}. (To detrend the data and perform aCompCor
#'  in the same regression, \code{nuisance} can be set to DCT bases obtained with
#'  the function \code{\link{dct_bases}}.)
#'
#' If \code{prelim_nuisance} includes an intercept, note that the PCs may retain
#'  a small component along the constant, because the centering of each noise
#'  ROI voxel is by the median (see \code{\link{noise_PCs}}).
#'
#' @inheritParams data_CompCor_Params
#' @inheritParams noise_Params
#' @param center,scale Center the columns of the noise ROI data by their medians,
#'  and scale by their MADs? Default: \code{TRUE} for both. Note that this argument
#'  affects the noise ROI data and not the data that is being cleaned with aCompCor.
#'  Centering and scaling of the data being cleaned can be done after this function call.
#' @param prelim_nuisance Nuisance signals to regress from the noise ROIs before
#'  computing the PCs. Should be a \eqn{T} by \eqn{P} numeric matrix where
#'  \eqn{P} represents the number of preliminary nuisance signals. Default:
#'  \code{NULL} (do not do any preliminary nuisance regression).
#' @param nuisance Nuisance signals to regress from each data column in addition
#'  to the noise ROI PCs. Should be a \eqn{T} by \eqn{N} numeric matrix where
#'  \eqn{N} represents the number of nuisance signals. To not perform any nuisance
#'  regression set this argument to \code{NULL}, \code{0}, or \code{FALSE}.
#'  Default: \code{NULL}.
#' @return A list with entries \code{"data"}, \code{"noise"}, and potentially
#'  \code{"ROI_data"}.
#'
#'  The entry \code{"data"} will be a matrix where each row or column
#'  corresponds to a data location (if it was originally an array, the locations
#'  will be voxels in spatial order). Each location's time series will have each
#'  noise PC regressed from it. This entry will be \code{NULL} if there was no
#'  data.
#'
#'  The entry \code{"noise"} is a \code{"CompCor_noise"} object: a list with the
#'  noise PC scores (\code{"PCs"}), their corresponding variance (\code{"var"}),
#'  and their ROI mask (\code{"ROI_noise"}), for each noise ROI.
#'
#'  If the data ROI is not all \code{TRUE}, the entry \code{"ROI_data"} will have
#'  the ROI mask for the data.
#'
#' @section References:
#'  \itemize{
#'    \item{Behzadi, Y., Restom, K., Liau, J. & Liu, T. T. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. NeuroImage 37, 90-101 (2007).}
#'    \item{Muschelli, J. et al. Reduction of motion-related artifacts in resting state fMRI using aCompCor. NeuroImage 96, 22-35 (2014).}
#' }
#'
#' @export
#' @seealso CompCor_clean, CompCor_HCP_calc
CompCor <- function(
  X, ROI_data="infer", ROI_noise=NULL,
  noise_nPC=5, noise_erosion=NULL,
  center=TRUE, scale=TRUE,
  prelim_nuisance=NULL, nuisance=NULL
  ){

  if (is.null(ROI_noise)) {stop("`ROI_noise` must be provided. (It is not inferred.)")}

  out1 <- format_data(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise,
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )

  # Compute the noise PCs (with the preliminary nuisance regression).
  out2 <- noise_PCs(
    X_noise=out1$X_noise, center=center, scale=scale,
    noise_nPC=out1$noise_nPC, prelim_nuisance=prelim_nuisance
  )
  noise <- new_CompCor_noise(out2$PCs, out2$var, out1$ROI_noise)

  # Clean the data with the noise PCs.
  data <- NULL
  if (any(out1$ROI_data)) {
    data <- CompCor_clean(out1$X, noise, nuisance=nuisance)
  }

  z <- list(data=data, noise=noise)
  if (!all(out1$ROI_data)) { z$ROI_data <- out1$ROI_data }
  class(z) <- "CompCor"
  z
}
