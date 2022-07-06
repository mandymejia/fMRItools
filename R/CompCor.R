#' CompCor: get noise components
#' 
#' Get noise components for aCompCor.
#'
#' @param X_noise The noise ROIs data
#' @param center,scale Center & scale robustly
#' @param noise_nPC Number of PCs to obtain for each noise ROI
#' 
#' @return A list with components X, X_noise, ROI_data, ROI_noise, noise_nPC,
#'  noise_erosion, noise_comps, and noise_var.
#' 
#' @keywords internal
CompCor.noise_comps <- function(X_noise, center, scale, noise_nPC){
  TOL <- 1e-8

  N <- length(X_noise)
  noise_comps <- vector("list", N); names(noise_comps) <- names(X_noise)
  noise_var <- vector("list", N); names(noise_var) <- names(X_noise)
  noise_vartotal <- vector("list", N); names(noise_vartotal) <- names(X_noise)

  if (length(noise_nPC) == 1) {
    noise_nPC <- as.list(rep(noise_nPC, N))
  } else {
    noise_nPC <- as.list(noise_nPC)
  }
  names(noise_nPC) <- names(X_noise)

  for (ii in 1:N) {
    T_ <- nrow(X_noise[[ii]])
    if (ncol(X_noise[[ii]])==0) { next }

    # Transpose.
    X_noise[[ii]] <- t(X_noise[[ii]])
    #	Center.
    if (center) { X_noise[[ii]] <- X_noise[[ii]] - c(rowMedians(X_noise[[ii]], na.rm=TRUE)) }
    # Compute MADs.
    mad <- 1.4826 * rowMedians(abs(X_noise[[ii]]), na.rm=TRUE)
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
    noise_comps[[ii]] <- x$u[,seq(noise_nPC[[ii]]),drop=FALSE]
    noise_var[[ii]] <- noise_var[[ii]][seq(noise_nPC[[ii]])]
  }

  list(noise_comps=noise_comps, noise_var=noise_var, noise_vartotal=noise_vartotal)
}

#' Anatomical CompCor
#'
#' The aCompCor algorithm for denoising fMRI data using noise ROIs data
#' 
#' First, the principal components (PCs) of each noise region of interest (ROI) 
#'  are calculated. For each ROI, voxels are centered and scaled 
#'  (can be disabled with the arguments \code{center} and \code{scale}), 
#'  and then the PCs are calculated via the singular value decomposition. 
#' 
#' Next, aCompCor is performed to remove the shared variation between the noise ROI
#'  PCs and each location in the data. This is accomplished by a nuisance regression
#'  using a design matrix with the noise ROI PCs, any additional regressors specified
#'  by \code{nuisance}, and an intercept term. (To detrend the data and perform aCompCor
#'  in the same regression, \code{nuisance} can be set to DCT bases obtained with
#'  the function \code{\link{dct_bases}}.)
#'
#' @inheritParams data_CompCor_Params
#' @inheritParams noise_Params
#' @param center,scale Center the columns of the noise ROI data by their medians, 
#'  and scale by their MADs? Default: \code{TRUE} for both. Note that this argument
#'  affects the noise ROI data and not the data that is being cleaned with aCompCor.
#'  Centering and scaling of the data being cleaned can be done after this function call.
#' @param nuisance Nuisance signals to regress from each data column in addition
#'  to the noise ROI PCs. Should be a \eqn{T} by \eqn{N} numeric matrix where 
#'  \eqn{N} represents the number of nuisance signals. To not perform any nuisance
#'  regression set this argument to \code{NULL}, \code{0}, or \code{FALSE}.
#'  Default: \code{NULL}.
#' @return A list with entries \code{"data"}, \code{"noise"}, and potentially
#'  \code{"ROI_data"}.
#'
#'  The entry \code{"data"} will be a \code{V x T} matrix where each row corresponds to a
#'  data location (if it was originally an array, the locations will be voxels in spatial
#'  order). Each row will be a time series with each noise PC regressed from it. This entry
#'  will be \code{NULL} if there was no data.
#'
#'  The entry \code{"noise"} is a list of noise PC scores, their corresponding variance,
#'  and their ROI mask, for each noise ROI.
#' 
#'  If the data ROI is not all \code{TRUE}, the entry \code{"ROI_data"} will have
#'  the ROI mask for the data.
#'
#' @importFrom robustbase rowMedians
#' 
#' @section References:
#'  \itemize{
#'    \item{Behzadi, Y., Restom, K., Liau, J. & Liu, T. T. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. NeuroImage 37, 90-101 (2007).}
#'    \item{Muschelli, J. et al. Reduction of motion-related artifacts in resting state fMRI using aCompCor. NeuroImage 96, 22-35 (2014).}
#' }
#'
#' @export
#' @seealso CompCor_HCP
CompCor <- function(
  X, ROI_data="infer", ROI_noise=NULL, 
  noise_nPC=5, noise_erosion=NULL, 
  center=TRUE, scale=TRUE, nuisance=NULL
  ){

  if (is.null(ROI_noise)) {stop("`ROI_noise` must be provided. (It is not inferred.)")}

  out1 <- format_data(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise, 
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )

  out2 <- CompCor.noise_comps(
    X_noise=out1$X_noise, center=center, scale=scale, noise_nPC=out1$noise_nPC
  )

  T_ <- nrow(out1$X)

  if (any(out1$ROI_data)) {
    # Perform nuisance regression.
    design <- do.call(cbind, out2$noise_comps)
    if (!(is.null(nuisance) || isFALSE(nuisance) || identical(nuisance, 0))) {
      if (nrow(nuisance) != T_) { stop("`nuisance` does not have the same number of timepoints as `X`.") }
      nuisance <- check_design_matrix(nuisance)
      design <- cbind(design, nuisance)
    }
    design <- check_design_matrix(cbind(1, design), T_)
    out1$X <- nuisance_regression(out1$X, design)

    # # Normalize data.
    # out1$X <- t(out1$X)
    # if (center) { out1$X <- out1$X - c(rowMedians(out1$X, na.rm=TRUE)) }
    # if (scale) { 
    #   mad <- 1.4826 * rowMedians(abs(out1$X), na.rm=TRUE)
    #   mad_inv <- ifelse(mad < 1e-8, 0, 1/mad)
    #   out1$X <- out1$X * mad_inv
    # }
    # out1$X <- t(out1$X)
  } else {
    out1["X"] <- list(NULL)
  }

  z <- list(
    data = out1$X, 
    noise = list(PCs=out2$noise_comps, var=out2$noise_var, ROI_noise=out1$ROI_noise)
  )
  if (!all(out1$ROI_data)) { z$ROI_data <- out1$ROI_data }
  class(z) <- "CompCor"
  z
}