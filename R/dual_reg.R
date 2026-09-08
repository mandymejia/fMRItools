#' Dual Regression
#'
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T}). Rows will be
#'  centered.
#' @param GICA Group-level independent components (\eqn{V \times Q})
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
# Below: inherit the scaling and temporal filtering parameters' documentation.
#' @inheritParams norm_BOLD 
#'
#' @return A list containing
#'  the subject-level independent components \strong{S} (\eqn{Q \times V}),
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#'
#' @export
#' @examples
#' nT <- 30
#' nV <- 400
#' nQ <- 7
#' mU <- matrix(rnorm(nV*nQ), nrow=nV)
#' mS <- mU %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ)
#' BOLD <- mS + rnorm(nV*nT, sd=.05) + 10
#' GICA <- mU
#' dual_reg(BOLD=BOLD, GICA=mU, scale_sm_FWHM=Inf, TR=.72)
#'
dual_reg <- function(
  BOLD, GICA,
  scale_by=c("mean", "sd", "none"),
  scale_sm_FWHM=4,
  scale_sm_xifti=NULL,
  TR=NULL, hpf=.01, lpf=NULL,
  GSR=FALSE){

  # [NOTE] to devs: if updating this function, please also make appropriate
  #   updates to `dual_reg_parc`.

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  scale_by <- match.arg(scale_by, c("mean", "sd", "none"))
  stopifnot(fMRItools::is_1(scale_sm_FWHM, "numeric"))
  # [NOTE]: 
  #   `scale_by=="none"` skips scaling completely
  #   `scale_sm=="none"` skips smoothing of scale estimates
  scale_sm <- switch(
    as.character(scale_sm_FWHM), 
    "0"="none", "Inf"="global", "local"
  )
  if (scale_sm == "local") {
    stopifnot(scale_sm_FWHM > 0)
    if (is.null(scale_sm_xifti)) {
      warning("Skipping smoothing of scale estimate because `scale_sm_xifti` ",
        "was not provided. If intended, set `scale_sm_FWHM=0` to disable this ",
        "warning.")
      scale_sm_FWHM <- 0; scale_sm <- "none"
    } else {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
      }
      stopifnot(ciftiTools::is.xifti(scale_sm_xifti))
    }
  }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (any(is.na(BOLD))) { stop("`NA` values in `BOLD` not supported with DR.") }
  if (any(is.na(GICA))) { stop("`NA` values in `GICA` not supported with DR.") }

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center each voxel timecourse. 
  #  Do not center the image at each timepoint unless `GSR == TRUE`.
  # Standardize scale if `scale_by != "none`, and do temporal filtering.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=GSR,
    scale_by=scale_by, scale_sm_FWHM=scale_sm_FWHM, 
    scale_sm_xifti=scale_sm_xifti, 
    # [NOTE]: could add the below arguments?
    # scale_sm_xifti_mask=scale_sm_xifti_mask, scale_precomp=scale_precomp,
    TR=TR, hpf=hpf, lpf=lpf
  ))

  # Center each group IC across space. (Used to be a function argument.)
  GICA <- colCenter(GICA)

  # Estimate A (IC timeseries).
  # We need to center `BOLD` across space because the linear model has no intercept.
  A <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% GICA) %*% chol2inv(chol(crossprod(GICA)))

  # Center each subject IC timecourse across time.
  # (Redundant. Since BOLD is column-centered, A is already column-centered.)
  # A <- colCenter(A)

  # Normalize each subject IC timecourse to constrain the ICA. (Used to be a function argument.)
  A <- scale(A)

  # Check rank of `A`.
  A_rank <- qr(A)$rank
  if (A_rank < ncol(A)) {
    warning(
      "DR has estimated an `A` matrix that has ", ncol(A), " columns, but its rank is ", A_rank, ". ",
      "An `A` matrix that is not full rank can occur when the number of group ICs approaches the number of volumes in the subject data. ",
      "This problem can be avoided by using a group ICA with fewer components, ",
      "or by providing more volumes of data. ",
      "Continuing, but an error may occur in further calculations."
    )
  }

  # Estimate S (IC maps).
  # Don't worry about the intercept: `BOLD` and `A` are centered across time.
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  # Re-estimate A (IC timeseries) based on the subject-level IC maps
  # We need to center `BOLD` across space because the linear model has no intercept.
  S_ctr <- colCenter(t(S))
  A2 <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% S_ctr) %*% chol2inv(chol(crossprod(S_ctr)))

  #return result
  list(S = S, A = A, A2 = A2)
}
