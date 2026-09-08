#' Normalize BOLD data
#'
#' Center the data across space and/or time, detrend, and scale, in that order.
#'  For dual regression, row centering is required and column centering is not
#'  recommended. Scaling and detrending depend on the user preference.
#'
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param center_rows,center_cols Center BOLD data across rows (each data
#'  location's time series) or columns (each time point's image)? Default:
#'  \code{TRUE} for row centering, and \code{FALSE} for column centering.
#' @param scale_by Scale the BOLD at each voxel based on either its 
#'  \code{"mean"} (default), its \code{"sd"}, or do not scale (\code{"none"}). 
#'  Mean scaling cannot be used if the \code{BOLD} have already been de-meaned. 
#' @param scale_sm_FWHM Full width at half maximum (FWHM) for smoothing the
#'  estimates of scale across brain locations (see \code{scale_by}), to reduce
#'  the variance of the estimates. Set to \code{0} to disable smoothing, or
#'  \code{Inf} for "global" smoothing (estimate and use one measure of scale 
#'  across the entire brain). Otherwise, for "local" smoothing, this should be a
#'  positive number. Default: \code{4}.
#' @param scale_sm_xifti Required only for "local" smoothing of the scale
#'  estimates with CIFTI data (see \code{scale_sm_FWHM}). To smooth the scale 
#'  estimates, provide a \code{"xifti"} object aligned with \code{BOLD}. If no
#'  \code{"xifti"} object is provided (default) smoothing must be skipped.
#' @param scale_sm_xifti_mask For local smoothing of scale estimates, the data
#'  must be unmasked to be mapped back to the surface. So if the data are 
#'  masked, provide the mask here.
#' @param scale_precomp (Optional) pre-computed image or scalar of BOLD average
#'  or SD, for \code{"mean"} or \code{"sd"} scaling respectively.
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' @param hpf,lpf The frequencies at which to apply temporal filtering to the
#'  data during pre-processing, in Hertz. Set either to \code{NULL} to disable.
#'  Default: \code{0.01} Hz highpass filter, and \code{NULL} for the lowpass 
#'  filter (disabled). Filtering is accomplished by nuisance regression of
#'  discrete cosine transform (DCT) bases.
#' 
#'  The highpass filter serves to detrend the data, since low-frequency 
#'  variance is associated with noise. The lowpass filter removes high-frequency
#'  variance, which is also thought to be from non-neuronal noise.
#' 
#'  Note the \code{TR} argument is required for temporal filtering. If
#'  \code{TR} is not provided, \code{hpf} and \code{lpf} will be ignored.
#'
#' @return Normalized BOLD data matrix (\eqn{V \times T})
#'
#' @export
#'
norm_BOLD <- function(
  BOLD, center_rows=TRUE, center_cols=FALSE,
  scale_by=c("mean", "sd", "none"),
  scale_sm_FWHM=4,
  scale_sm_xifti=NULL,
  scale_sm_xifti_mask=NULL, 
  scale_precomp=NULL,
  TR=NULL, hpf=.01, lpf=NULL){

  nT <- ncol(BOLD)
  nV <- nrow(BOLD)
  if (nT > nV) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  scale_by <- match.arg(scale_by, c("mean", "sd", "none"))
  stopifnot(fMRItools::is_1(scale_sm_FWHM, "numeric"))
  # [NOTE]: 
  #   `scale_by=="none"` skips scaling completely
  #   `scale_sm=="none"` skips smoothing of scale estimates
  scale_sm <- switch(
    as.character(scale_sm_FWHM), 
    "0"="none", "Inf"="global", "local"
  )
  if (scale_by != "none" && scale_sm == "local") {
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
      if (!is.null(scale_sm_xifti_mask)) {
        stopifnot(is.vector(scale_sm_xifti_mask) && is.logical(scale_sm_xifti_mask))
        stopifnot(sum(scale_sm_xifti_mask) == nV)
      }
    }
  }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (length(hpf)==1 && hpf==0) { hpf <- NULL }
  if (length(lpf)==1 && lpf==Inf) { lpf <- NULL }
  if (is.null(TR)) {
    if (!is.null(hpf)) {
      if (hpf==.01) {
        message("Setting `hpf=NULL` because `TR` was not provided. Either provide `TR` or set `hpf=NULL` to disable this message.")
        hpf <- NULL
      } else {
        stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=NULL`.")
      }
    }
    if (!is.null(lpf)) {
      stop("Cannot apply `lpf` because `TR` was not provided. Either provide `TR` or set `lpf=NULL`.")
    }
  } else {
    stopifnot(is_posNum(TR))
    stopifnot(is.null(hpf) || is_posNum(hpf, zero_ok=TRUE))
    stopifnot(is.null(lpf) || is_posNum(lpf))
  }

  # Get `voxMeans` before doing anything else. ---------------------------------
  voxMeans <- rowMeans(BOLD, na.rm=TRUE)

  # Center. --------------------------------------------------------------------
  if (center_rows || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD - voxMeans)
    } else {
      BOLD <- t(BOLD)
    }
    # Center each image (across space).
    if (center_cols) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      BOLD <- t(BOLD)
    }
  }

  # Apply the temporal filter. -------------------------------------------------
  # [TO DO]: consider using `fsl_bptf` for LPF rather than DCT (too many bases?)
  # [NOTE]: If `center_cols`, columns won't be exactly centered anymore after the filter.
  if (!is.null(hpf) || !is.null(lpf)) {
    dct <- fMRItools::temporal_filter(
      X=ncol(BOLD), TR=TR, hpf=hpf, lpf=lpf, method="DCT", verbose=FALSE
    ) # [TO DO] carry over verbose arg?
    if (!is.null(dct) && nrow(dct) > 0) {
      BOLD <- nuisance_regression(BOLD, cbind(1, dct))
      if (!center_rows) { BOLD <- BOLD + voxMeans }
    }
  }

  # Scale. ---------------------------------------------------------------------
  if (scale_by == "none") { return(invisible(BOLD)) }

  ## Get scale measure.
  if (!is.null(scale_precomp)) {
    if (scale_sm %in% c("none", "local")) {
      stopifnot(length(scale_precomp) == nV)
    }
    scale_meas <- scale_precomp
  } else {
    scale_meas <- switch(scale_by,
      mean = voxMeans,
      sd = sqrt(rowVars(BOLD, na.rm=TRUE))
    )
  }

  if (mean(scale_meas, na.rm=TRUE) < 1e-8) {
    stop("Estimated scale is near zero.")
  }

  ## Smooth scale measure, if applicable.
  ### Global.
  if (scale_sm == "global") { scale_meas <- mean(scale_meas, na.rm=TRUE) }

  ### Local.
  if (scale_sm == "local") {
    # Check `scale_sm_xifti` is valid.
    is_masked <- !is.null(scale_sm_xifti_mask)

    # Un-mask, if applicable.
    if (is_masked) {
      scale_meas <- c(unmask_mat(as.matrix(scale_meas), scale_sm_xifti_mask))
      nV <- length(scale_meas)
    }
    if (nV != nrow(scale_sm_xifti)) {
      stop("`scale_sm_xifti` not compatible with `BOLD`: different spatial dimensions.")
    }
    if (!is.null(scale_sm_xifti$meta$cifti$intent) && scale_sm_xifti$meta$cifti$intent == 3007) {
      scale_sm_xifti <- ciftiTools::convert_xifti(scale_sm_xifti, "dscalar")
    }

    # Convert `scale_meas` to `"xifti"`.
    scale_meas <- ciftiTools::newdata_xifti(ciftiTools::select_xifti(scale_sm_xifti, 1), scale_meas)
    scale_meas <- ciftiTools::move_to_mwall(scale_meas, NA)
    if (!is.null(scale_meas$data$subcort)) {
      sub_mask <- !is.na(scale_meas$data$subcort[,1])
      scale_meas$data$subcort <- scale_meas$data$subcort[sub_mask,,drop=FALSE]
      scale_meas$meta$subcort$labels <- scale_meas$meta$subcort$labels[sub_mask]
      scale_meas$meta$subcort$mask[scale_meas$meta$subcort$mask][!sub_mask] <- FALSE
    }

    # Smooth `scale_meas`.
    scale_meas <- ciftiTools::smooth_xifti(scale_meas, surf_FWHM=scale_sm_FWHM, vol_FWHM=scale_sm_FWHM)
    
    # Convert `scale_meas` back to matrix.
    scale_meas <- c(as.matrix(scale_meas))

    # Re-mask, if applicable, to match `BOLD`'s original (masked) spatial dimension.
    if (is_masked) {
      scale_meas <- scale_meas[scale_sm_xifti_mask]
    }
  }

  # Apply scaling.
  BOLD <- BOLD / scale_meas

  if (scale_by == "mean") {
    # Multiply by 100 to make it "percent change of signal".
    BOLD <- BOLD * 100
  }

  invisible(BOLD)
}
