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
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation.
#' @param scale_sm_xifti,scale_sm_FWHM Only applies if \code{scale=="local"} and
#'  \code{BOLD} represents CIFTI-format data. To smooth the standard deviation
#'  estimates used for local scaling, provide a \code{"xifti"} object with data
#'  locations in alignment with \code{BOLD}, as well as the smoothing FWHM
#'  (default: \code{2}). If no \code{"xifti"} object is provided (default), do
#'  not smooth.
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' @param hpf The frequency at which to apply a highpass filter to the data
#'  during pre-processing, in Hertz. Default: \code{0.01} Hertz. Set to \code{0}
#'  to disable the highpass filter.
#' 
#' 
#'  The highpass filter serves to detrend the data, since low-frequency 
#'  variance is associated with noise. Highpass filtering is accomplished by 
#'  nuisance regression of discrete cosine transform (DCT) bases. 
#' 
#'  Note the \code{TR} argument is required for highpass filtering. If
#'  \code{TR} is not provided, \code{hpf} will be ignored.
#'
#' @return Normalized BOLD data matrix (\eqn{V \times T})
#'
#' @importFrom matrixStats rowVars
#' @export
#'
norm_BOLD <- function(
  BOLD, center_rows=TRUE, center_cols=FALSE,
  scale=c("local", "global", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01){

  nT <- ncol(BOLD)
  nV <- nrow(BOLD)
  if (nT > nV) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  if (!is.null(scale_sm_xifti)) {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
    stopifnot(ciftiTools::is.xifti(scale_sm_xifti))
  }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (is.null(hpf)) { hpf <- 0 }
  if (is.null(TR)) {
    if (hpf==.01) {
      message("Setting `hpf=0` because `TR` was not provided. Either provide `TR` or set `hpf=0` to disable this message.")
      hpf <- 0
    } else if (hpf!=0) {
      stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=0`.")
    }
  } else {
    stopifnot(is_posNum(TR))
    stopifnot(is_posNum(hpf, zero_ok=TRUE))
  }

  # Center.
  if (center_rows || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
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

  # Detrend.
  # [NOTE]: If `center_cols`, columns won't be centered anymore after detrending.
  if (hpf > 0) {
    nDCT <- round(dct_convert(T_=nT, TR=TR, f=hpf))
    if (nDCT == 0) {
      warning("For the low `hpf` and at the data TR and length, the closest number of DCT bases to use for detrending is zero. Using one instead. See `dct_convert`.")
      nDCT <- 1
    }
    if (!center_rows) { voxMeans <- rowMeans(BOLD, na.rm=TRUE) }
    BOLD <- nuisance_regression(
      BOLD,
      cbind(1, dct_bases(nT, nDCT))
    )
    if (!center_rows) { BOLD <- BOLD + voxMeans }
  }

  # Scale.
  # Get scale at each location.
  if (scale != "none") { sig <- sqrt(matrixStats::rowVars(BOLD, na.rm=TRUE)) }
  # Global scaling: take mean scale across all locations, and use that.
  if (scale == "global") {
    sig <- mean(sig, na.rm=TRUE)
    if (sig < 1e-8) {
      warning("Estimated scale is near zero. Skipping scaling.")
    } else {
      # Apply global scaling.
      BOLD <- BOLD / sig
    }
  # Local scaling: use estimate of scale at each location.
  } else if (scale == "local") {
    # Smooth estimates, if applicable.
    if (!is.null(scale_sm_xifti) && (scale_sm_FWHM != 0)) {
      # Check `scale_sm_xifti` is valid.
      if (nV != nrow(scale_sm_xifti)) {
        stop("`scale_sm_xifti` not compatible with `BOLD`: different spatial dimensions.")
      }
      if (!is.null(scale_sm_xifti$meta$cifti$intent) && scale_sm_xifti$meta$cifti$intent == 3007) {
        scale_sm_xifti <- ciftiTools::convert_xifti(scale_sm_xifti, "dscalar")
      }
      # Compute and smooth the SD.
      sig <- ciftiTools::newdata_xifti(ciftiTools::select_xifti(scale_sm_xifti, 1), sig)
      sig <- ciftiTools::smooth_xifti(sig, surf_FWHM=scale_sm_FWHM, vol_FWHM=scale_sm_FWHM)
      sig <- c(as.matrix(sig))
    }
    # Apply local scaling.
    BOLD <- BOLD / sig
  }

  BOLD
}
