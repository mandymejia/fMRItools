#' Normalize BOLD data
#'
#' Clean the data (scrubbing, detrending, nuisance regression, and temporal
#'  filtering in a simultaneous framework), then center it, then scale it, in
#'  that order.
#'
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param drop_first (Optional) Number of volumes to drop from the start of each
#'  BOLD session. Default: \code{0}.
#' @param nuisance (Optional) Nuisance matrix to regress from the BOLD data.
#'
#'  Nuisance regression is performed in a simultaneous regression with any spike
#'  regressors from \code{scrub} and DCT bases from \code{hpf}.
#'
#'  Note that the nuisance matrix should be provided with timepoints matching
#'  the original \code{BOLD} irregardless of \code{drop_first}.
#'  Nuisance matrices will be truncated automatically if \code{drop_first>0}.
#' @param scrub (Optional) Numeric vector of integers giving the indices
#'  of volumes to scrub from the BOLD data. (List the volumes to remove, not the
#'  ones to keep.)
#'
#'  Scrubbing is performed within a nuisance regression by adding a spike
#'  regressor to the nuisance design matrix for each volume to scrub.
#'
#'  Note that indices are counted beginning with the first index in the
#'  \code{BOLD} session irregardless of \code{drop_first}. The indices will be
#'  adjusted automatically if \code{drop_first>0}.
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' @param hpf,lpf The frequencies at which to apply temporal filtering to the
#'  data during pre-processing, in Hertz. Set either to \code{NULL} to disable.
#'  Default: \code{0.01} Hz highpass filter, and \code{NULL} for the lowpass
#'  filter (disabled). Filtering is accomplished by nuisance regression of
#'  discrete cosine transform (DCT) bases.
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
#' @param scale_FUN If \code{scale_by} is \code{"FUN"}, this function is used
#'  to compute the scaling measure from BOLD, after it has been centered
#'  (and so after dropping any volumes, and applying nuisance regression,
#'  scrubbing, and temporal filtering). Otherwise, ignored.
#' @param give_stats Return the intercept and residual SD estimates from the
#'  nuisance regression? Default: \code{FALSE}
#'
#'  The highpass filter serves to detrend the data, since low-frequency
#'  variance is associated with noise. The lowpass filter removes high-frequency
#'  variance, which is also thought to be from non-neuronal noise.
#'
#'  Note the \code{TR} argument is required for temporal filtering. If
#'  \code{TR} is not provided, \code{hpf} and \code{lpf} will be ignored.
#'
#' @return Normalized BOLD data matrix (\eqn{V \times T}), or if \code{give_stats},
#'  a list with three elements: the normed BOLD, the intercept estimate, and the
#'  residual SD estimate.
#'
#' @export
#'
norm_BOLD <- function(
  BOLD, drop_first=0,
  nuisance=NULL, scrub=NULL,
  TR=NULL, hpf=.01, lpf=NULL,
  center_rows=TRUE, center_cols=FALSE,
  scale_by=c("mean", "sd", "FUN", "none"),
  scale_sm_FWHM=4,
  scale_sm_xifti=NULL, scale_sm_xifti_mask=NULL,
  scale_FUN=NULL,
  give_stats=FALSE){

  # Check arguments. -----------------------------------------------------------
  nV <- nrow(BOLD)
  nT <- ncol(BOLD)
  if (nV < nT) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)

  scale_by <- match.arg(scale_by, c("mean", "sd", "FUN", "none"))
  stopifnot(is_1(scale_sm_FWHM, "numeric"))
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
      if (dim(scale_sm_xifti)[1] == 0) {
        if (!is.null(scale_sm_xifti$surf$cortex_left)) {
          scale_sm_xifti$data$cortex_left <- as.matrix(rep(0, nrow(scale_sm_xifti$surf$cortex_left$vertices)))
        }
        if (!is.null(scale_sm_xifti$surf$cortex_right)) {
          scale_sm_xifti$data$cortex_right <- as.matrix(rep(0, nrow(scale_sm_xifti$surf$cortex_right$vertices)))
        }
      }
      if (!is.null(scale_sm_xifti_mask)) {
        stopifnot(is.vector(scale_sm_xifti_mask) && is.logical(scale_sm_xifti_mask))
        stopifnot(sum(scale_sm_xifti_mask) == nV)
      }
    }
  }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)

  # Ensure `nuisance` is a numeric matrix.
  if (!is.null(nuisance)) {
    stopifnot(is.numeric(nuisance) && is.matrix(nuisance) && nrow(nuisance)==nT)
    if (ncol(nuisance) == 0) { nuisance <- NULL }
  }

  # Create `scrub_mat` (spike regressor matrix) if any `scrub`.
  # Do not scrub volumes that are one of the first `drop_first`.
  if (!is.null(scrub)) {
    stopifnot(is.numeric(scrub))
    if (!is.null(drop_first)) {
      stopifnot(is_posNum(drop_first, zero_ok=TRUE))
      scrub <- scrub[scrub > drop_first] - drop_first
    }
  }
  if (length(scrub) == 0) { scrub <- NULL }

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

  stopifnot(is_1(give_stats, "logical"))

  # `drop_first`. --------------------------------------------------------------
  stopifnot(is_posNum(drop_first, zero_ok=TRUE))
  if (drop_first > 0) {
    stopifnot(drop_first < nT - 2)
    # Drop columns from BOLD; drop rows from nuisance.
    # (Already done for scrubbing.)
    BOLD <- BOLD[,-seq(drop_first),drop=FALSE]
    if (!is.null(nuisance)) {
      nuisance <- nuisance[-seq(drop_first),,drop=FALSE]
    }
    # Recalculate nT.
    nT <- ncol(BOLD)
  }

  # Nuisance regression. -------------------------------------------------------
  # Includes any: input regressors; scrubbing; temporal filtering.
  # Note: HPF is done using DCT bases in the regression, but
  #       LPF is done using `fsl_bptf` applied to the data + nuisance matrices
  #       prior the regression. (For a linear filter, that would be equivalent)
  #       to simultaneous regression.)
  # Also, calculate the intercept and residuals SD for returning.

  ### Prepare the BOLD and design matrix. --------------------------------------
  add_to_nuis <- function(x, nuis) {
    if (is.null(nuis)) {
      if (is.matrix(x)) { x } else { as.matrix(x, ncol=1) }
    } else {
      # Add new columns on the right
      cbind(nuis, x)
    }
  }

  # Intercept column
  big_nmat <- add_to_nuis(rep(1, nT), NULL)

  # HPF
  if (!is.null(hpf)) {
    nDCT <- round(dct_convert(nT, TR=TR, f=hpf))
    if (nDCT > 0) {
      big_nmat <- add_to_nuis(dct_bases(nT, nDCT), big_nmat)
    }
  }

  # `nuisance` (input regressors)
  if (!is.null(nuisance)) {
    stopifnot(is.numeric(nuisance) && is.matrix(nuisance) && nrow(nuisance)==nT)
    if (ncol(nuisance) > 0) {
      big_nmat <- add_to_nuis(nuisance, big_nmat)
    }
  }

  # Scrubbing
  if (!is.null(scrub)) {
    scrub_mat <- flags2spikes(scrub, nT)
    big_nmat <- add_to_nuis(scrub_mat, big_nmat)
  }

  ### LPF both `BOLD` and `big_nmat.` ------------------------------------------
  if (!is.null(lpf)) {
    LP_sigma <- 1 / (18 * lpf * TR)
    BOLD <- t(fsl_bptf(t(BOLD), HP_sigma=NULL, LP_sigma=LP_sigma))
    big_nmat <- fsl_bptf(big_nmat, HP_sigma=NULL, LP_sigma=LP_sigma)
  }

  ### Do the regression. -------------------------------------------------------
  # Calculate mu as the intercept from nuisance regression.
  BOLD_mu <- (solve(crossprod(big_nmat)) %*% t(big_nmat) %*% t(BOLD))[1,]

  # Do the nuisance regression. (BOLD is now centered.)
  BOLD <- nuisance_regression(BOLD, big_nmat)

  # Compute the DOF for SD estimate calculation.
  # (DOF lost by scrubbing accounted for in `nmat` rather than dropped columns).
  BOLD_dof <- nT - qr(big_nmat)$rank

  # Compute the SD estimates for SD scaling.
  BOLD_sd <- sqrt(rowSums(BOLD^2, na.rm=TRUE) / BOLD_dof)

  # Drop scrubbed volumes.
  if (!is.null(scrub)) { BOLD <- BOLD[,-scrub,drop=FALSE] }

  # Define the scaling measure.
  scale_meas <- switch(scale_by,
    "mean"=BOLD_mu,
    "sd"=BOLD_sd,
    "FUN"=NULL,
    "none"=NULL
  )

  # Center. --------------------------------------------------------------------
  if ((!center_rows) || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD)
    } else {
      BOLD <- t(BOLD + BOLD_mu)
    }
    # Center each image (across space).
    if (center_cols) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      BOLD <- t(BOLD)
    }
  }

  # Scale. ---------------------------------------------------------------------

  # Skip and return if no scaling.
  if (scale_by == "none") {
    out <- if (give_stats) {
      list(BOLD=BOLD, mu=BOLD_mu, sd=BOLD_sd)
    } else {
      BOLD
    }
    return(invisible(out))
  }

  if (scale_by == "FUN") {
    stopifnot(is.function(scale_FUN))
    scale_meas <- scale_FUN(BOLD)
  }

  if (mean(scale_meas, na.rm=TRUE) < 1e-8) {
    stop("Estimated mean scale is near zero. ",
      "Set `scale = 'none'` or provide non-centered data.")
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
    scale_meas <- ciftiTools::convert_to_dscalar(scale_meas)
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

  # Checks.
  if (any(scale_meas < 0) || any(abs(scale_meas) < 1e-8)) {
    stop("Some locations have zero or negative scaling measures. ",
      "Set `scale = 'none'` or double-check the data. Note, ",
      "mean scaling requires uncentered data, and ",
      "SD scaling requires that constant volumes be masked prior to `norm_BOLD`."
    )
  }

  if (scale_by == "mean") {
    # Effectively multiplies BOLD by 100 to make it "percent change of signal".
    scale_meas <- scale_meas / 100
  }

  # Apply scaling.
  BOLD <- BOLD / scale_meas

  # Return. --------------------------------------------------------------------
  out <- if (give_stats) {
    list(BOLD=BOLD, mu=BOLD_mu, sd=BOLD_sd, scale_meas=scale_meas)
  } else {
    BOLD
  }
  invisible(out)
}
