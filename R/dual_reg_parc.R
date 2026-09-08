#' Multiple regression for parcel data
#'
#' @param parc The parcellation as an integer vector.
#' @param parc_vals The parcel values (keys) in desired order, e.g.
#'  \code{sort(unique(parc))}.
#' @inheritParams norm_BOLD
#' @inheritParams dual_reg
#'
#' @return A list containing
#'  the subject-level independent components \strong{S} (\eqn{Q \times V}),
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#'
#' @importFrom matrixStats rowMedians
#' @export
#'
dual_reg_parc <- function(
  BOLD, parc, parc_vals,
  scale_by=c("mean", "sd", "none"),
  scale_sm_FWHM=4,
  scale_sm_xifti=NULL,
  TR=NULL, hpf=.01, lpf=NULL,
  GSR=FALSE){

  # [NOTE] to devs: if updating this function, please also make appropriate
  #   updates to `dual_reg`.

  stopifnot(is.matrix(BOLD))
  stopifnot(is.numeric(parc))
  parc <- as.matrix(parc)
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
  if (any(is.na(parc))) { stop("`NA` values in `parc` not supported with DR.") }

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(parc)) {
    stop('The number of voxels in dat (', nV, ') and parc (', nrow(parc), ') must match')
  }

  stopifnot(all(unique(c(parc)) %in% parc_vals))
  stopifnot(all(parc_vals %in% parc))
  nQ <- length(parc_vals) #number of parcels
  if(nQ > nV) warning('More parcels than voxels. Are you sure?')
  if(nQ > nT) warning('More parcels than time points. Are you sure?')

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

  # Estimate A (parcel timeseries).
  # Neccesary to temporarily center BOLD (like in standard dual regression since there is no intecept)
  BOLD_rowMeans <- rowMeans(BOLD, na.rm=TRUE)
  A <- matrix(NA, nrow=nT, ncol=nQ)
  for (qq in seq(nQ)) {
    BOLD_qq <- BOLD[,parc==parc_vals[qq],drop=FALSE]
    A[,qq] <- matrixStats::rowMedians(BOLD_qq - BOLD_rowMeans)
  }
  rm(BOLD_rowMeans)

  # Normalize each subject parcel timecourse. (Used to be a function argument.)
  normA <- TRUE
  if (normA) { A <- scale(A) }

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
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  # Re-estimate A (IC timeseries) based on the subject-level IC maps
  # We need to center `BOLD` across space because the linear model has no intercept.
  S_ctr <- colCenter(t(S))
  A2 <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% S_ctr) %*% chol2inv(chol(crossprod(S_ctr)))

  #return result
  list(S = S, A = A, A2 = A2)
}
