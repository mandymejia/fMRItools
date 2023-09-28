#' ICA-based harmonization
#'
#' Harmonize data using ICA.
#'
#' Let \eqn{V} be the number of data locations; \eqn{Q} be the number of group
#'  ICs, and \eqn{N} be the number of fMRI sessions.
#'
#
#  Note for developers: this function is modeled after \code{estimate_template}
#  in \code{templateICAr}.
#
#' @param BOLD Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, GIFTI file paths,
#'  \code{"gifti"} objects, NIFTI file paths, \code{"nifti"} objects,
#'  or \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data
#'  locations and \eqn{T} is the number of timepoints.
#'  @param ghemi Which hemisphere, in the case of a single GIFTI file?
#
#   If GIFTI or \code{"gifti"}, the input can also be a length two list,
#   where the first list element is a length \eqn{N} vector for the left
#   hemisphere and the second list element is a length \eqn{N} vector for the
#   right hemisphere.

#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also
#'  be a (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of
#'  \code{BOLD}. Its columns will be centered.
#' @param mask Required if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects, optional for other formats. For NIFTI, this is a brain map formatted as a
#'  binary array of the same spatial dimensions as the fMRI data, with
#'  \code{TRUE} corresponding to in-mask voxels. For other formats, a logical vector.
# @param inds Numeric indices of the group ICs to include in the template. If
#  \code{NULL}, use all group ICs (default).
#
#  If \code{inds} is provided, the ICs not included will be removed after
#  calculating dual regression, not before. This is because removing the ICs
#  prior to dual regression would leave unmodelled signals in the data, which
#  could bias the templates.
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation.
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents surface data (CIFTI or
#'  GIFTI). To smooth the standard deviation estimates used for local scaling,
#'  provide the surface geometries along which to smooth as GIFTI geometry files
#'  or \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
#'
#'  If \code{scale_sm_FWHM==0}, no smoothing of the local standard deviation
#'  estimates will be performed.
#'
#'  If \code{scale_sm_FWHM>0} but \code{scale_sm_surfL} and
#'  \code{scale_sm_surfR} are not provided, the default inflated surfaces from
#'  the HCP will be used.
#'
#'  To create a \code{"surf"} object from data, see
#'  \code{\link[ciftiTools]{make_surf}}. The surfaces must be in the same
#'  resolution as the \code{BOLD} data.
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
#' @param center_Bcols Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' @param maskTol For computing the dual regression results for each subject:
#'  tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  a subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for both scans and then
#'  the intersection of the masks is used, for each subject.
#' @param missingTol For computing the variance decomposition across all subjects:
#'  tolerance for number of subjects masked out due to low variance or missing
#'  values at a given location. If more than this many subjects are masked out,
#'  the location's value will be \code{NA} in the templates. \code{missingTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of subjects can be masked out
#'  at a given location.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' @keywords internal
#'
harmonize <- function(
  BOLD, ghemi=NULL,
  GICA, mask=NULL, #inds=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  center_Bcols=FALSE,
  brainstructures=c("all"),
  varTol=1e-6, maskTol=.1, missingTol=.1,
  verbose=TRUE){

  if (!requireNamespace("templateICAr", quietly = TRUE)) {
    stop("Package \"templateICAr\" needed. Please install it", call. = FALSE)
  }

  # Check arguments ------------------------------------------------------------

  # Simple argument checks.
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='local'`. Use `'local'` or `'global'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "local"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  stopifnot(fMRItools:::is_1(scale_sm_FWHM, "numeric"))
  if (is.null(hpf)) { hpf <- 0 }
  if (is.null(TR)) {
    if (hpf==.01) {
      message("Setting `hpf=0` because `TR` was not provided. Either provide `TR` or set `hpf=0` to disable this message.")
      hpf <- 0
    } else {
      stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=0`.")
    }
  } else {
    stopifnot(fMRItools:::is_posNum(TR))
    stopifnot(fMRItools:::is_posNum(hpf, zero_ok=TRUE))
  }
  stopifnot(fMRItools:::is_1(center_Bcols, "logical"))
  stopifnot(fMRItools:::is_1(varTol, "numeric"))
  if (varTol < 0) { cat("Setting `varTol=0`."); varTol <- 0 }
  stopifnot(fMRItools:::is_posNum(maskTol))
  stopifnot(fMRItools:::is_posNum(missingTol))
  stopifnot(fMRItools:::is_1(verbose, "logical"))

  # `BOLD` format --------------------------------------------------------------
  format <- fMRItools:::infer_format_ifti_vec(BOLD)[1]
  FORMAT <- fMRItools:::get_FORMAT(format)
  FORMAT_extn <- switch(FORMAT,
    CIFTI=".dscalar.nii",
    GIFTI=".func.gii",
    NIFTI=".nii",
    MATRIX=".rds"
  )
  nN <- length(BOLD)

  fMRItools:::check_req_ifti_pkg(FORMAT)

  # If BOLD is a CIFTI, GIFTI, NIFTI, or RDS file, check that the file paths exist.
  if (format %in% c("CIFTI", "GIFTI", "NIFTI", "RDS")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (any(missing_BOLD)) {
      warning(
        'There are ', missing_BOLD,
        ' scans in `BOLD` that do not exist. ',
        'These scans will be excluded from template estimation.'
      )
      BOLD <- BOLD[!missing_BOLD]
    }
  }

  # Check `scale_sm_FWHM`
  if (scale_sm_FWHM !=0 && FORMAT %in% c("NIFTI", "MATRIX")) {
    if (scale_sm_FWHM==2) {
      cat("Setting `scale_sm_FWHM == 0`.\n")
    } else {
      if (FORMAT == "NIFTI") {
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for volumetric data.).\n")
      } else {
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for data matrices: use CIFTI/GIFTI files.).\n")
      }
    }
    scale_sm_FWHM <- 0
  }

  # `GICA` ---------------------------------------------------------------------
  # Convert `GICA` to a numeric data matrix or array.
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) { GICA <- ciftiTools::read_cifti(GICA, brainstructures=brainstructures) }
    if (ciftiTools::is.xifti(GICA, messages=FALSE)) {
      xii1 <- ciftiTools::select_xifti(GICA, 1) # for formatting output
      GICA <- as.matrix(GICA)
    } else {
      # Get `xii1` from first data entry.
      xii1 <- BOLD[[1]]
      if (is.character(xii1)) {
        xii1 <- ciftiTools::read_cifti(xii1, brainstructures=brainstructures, idx=1)
      }
      xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(xii1, 1), "dscalar")
    }
    stopifnot(is.matrix(GICA))
  } else if (FORMAT == "GIFTI") {
    if (is.character(GICA)) { GICA <- gifti::readgii(GICA) }
    ghemi <- GICA$file_meta["AnatomicalStructurePrimary"]
    if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    GICA <- do.call(cbind, GICA$data)
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- RNifti::readNifti(GICA) }
    stopifnot(length(dim(GICA)) > 1)
  } else if (FORMAT == "MATRIX") {
    if (is.character(GICA)) { GICA <- readRDS(GICA) }
    stopifnot(is.matrix(GICA))
  }
  nQ <- dim(GICA)[length(dim(GICA))]

  # # `inds`.
  # if (!is.null(inds)) {
  #   if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
  #   nL <- length(inds)
  # } else {
  #   inds <- seq(nQ)
  #   nL <- nQ
  # }

  # [TO DO]: NA in GICA?

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array.
  # Check `GICA` and `mask` dimensions match.
  # Append NIFTI header from GICA to `mask`.
  # Vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- RNifti::readNifti(mask); mask <- array(as.logical(mask), dim=dim(mask)) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) {
      cat("Coercing `mask` to a logical array.\n")
      mask <- array(as.logical(mask), dim=dim(mask))
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[seq(length(dim(GICA))-1)] == nI))
      }
      # Append NIFTI header.
      mask <- RNifti::asNifti(array(mask, dim=c(dim(mask), 1)), reference=GICA)
      # Vectorize `GICA`.
      if (all(dim(GICA)[seq(length(dim(GICA))-1)] == nI)) {
        GICA <- matrix(GICA[rep(as.logical(mask), nQ)], ncol=nQ)
        stopifnot(nrow(GICA) == nV)
      }
    }
  } else {
    if (!is.null(mask)) {
      #warning("Ignoring `mask`, which is only applicable to NIFTI data.")
      #mask <- NULL
      nI <- length(mask); nV <- sum(mask)
    } else {
      nI <- nV <- nrow(GICA)
    }
  }

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- fMRItools:::colCenter(GICA) }

  # Print summary of data ------------------------------------------------------
  format2 <- if (format == "data") { "numeric matrix" } else { format }
  if (verbose) {
    cat("Data input format:    ", format2, "\n")
    cat("Image dimensions:     ", paste(nI, collapse=" x "), "\n")
    cat('Masked locations:     ', nV, "\n")
    if (FORMAT == "GIFTI") {
      cat("Cortex hemisphere:    ", ghemi, "\n")
    }
    cat('Number of group ICs:  ', nQ, "\n")
    cat('Number of sessions:   ', nN, "\n")
  }

  # Dual regression ------------------------------------------------------------

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- fMRItools:::colCenter(GICA) }

  S0 <- array(0, dim = c(nN, nQ, nV))
  A0 <- vector("list", nN)
  G0 <- array(NA, dim = c(nN, nQ, nQ))
  for (ii in seq(nN)) {
    if (verbose) { cat(paste0(
      '\nReading and analyzing data for subject ', ii,' of ', nN, '.\n'
    )) }

    DR0_ii <- harmonize_DR_oneBOLD(
      BOLD[[ii]], mask=mask,
      ghemi=ghemi,
      format=format,
      GICA=GICA,
      center_Bcols=center_Bcols,
      scale=scale,
      scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
      scale_sm_FWHM=scale_sm_FWHM,
      TR=TR, hpf=hpf,
      brainstructures=brainstructures,
      varTol=varTol, maskTol=maskTol,
      verbose=verbose
    )
    S0[ii,,] <- DR0_ii$S
    A0[[ii]] <- DR0_ii$A
    G0[ii,,] <- cov(DR0_ii$A)
  }

  # Process the features -------------------------------------------------------

  ### Use SVD to reduce dimensions of S_q <- MAYBE NOT!  USE SPATIAL MODEL INSTEAD?

  U0 <- V0 <- vector('list', nQ)
  for(qq in 1:nQ){

    if (verbose) { cat(paste0(
      '\nPerforming PCA on IC ', qq,' of ', nQ, '.\n'
    )) }

    #do PCA separately for each IC q
    S0_q <- S0[,qq,] #NxV

    #center across sessions so X'X is covariance matrix (add back later)
    #S0_q_mean <- colMeans(S0_q)
    S0_q <- t(S0_q) #- S0_q_mean #vector will be recycled by column

    #want to obtain U (NxP) where P << V
    SSt_q <- crossprod(S0_q) # S0 is currently VxN, so S0'S0 is NxN.
    svd_q <- svd(SSt_q, nv=0) #want to get U from SVD of S0' = UDV'. SVD of S0'S0 = UDV'VDU' = U D^2 U'

    #keep components explaining 99% of variation
    nP <- min(which(cumsum(svd_q$d)/sum(svd_q$d) > 0.99))
    U0[[qq]] <- svd_q$u[,1:nP]

    #visualize V' = (1/D)U'S0' to make sure they are sensible
    V0[[qq]] <- diag(1/svd_q$d[1:nP]) %*% t(svd_q$u[,1:nP]) %*% t(S0_q)
  }
  save(U0, file=file.path(dir_features, 'U.RData'))
  save(V0, file=file.path(dir_features, 'V.RData'))
  U0_1 <- U0[[1]]
  save(U0, file=file.path(dir_features, 'U.RData'))
  save(U0_1, file=file.path(dir_features, 'U1.RData'))
  save(S0, file=file.path(dir_features, 'S.RData'))
  S0_1 <- S0[,1,]
  save(S0_1, file=file.path(dir_features, 'S1.RData'))
  save(A0, file=file.path(dir_features, 'A.RData'))
  save(G0, file=file.path(dir_features, 'G.RData'))

  ### Project elements of G to a tangent space

  if (verbose) { cat('\nProjecting covariance matrices to tangent space.\n' }

  # Calculate the element-wise average of the covariance matrices
  G_avg <- apply(G0, c(2, 3), mean)

  # Project each covariance matrix to tangent space at base point defined by the Euclidean mean
  system.time(G_tangent <- apply(G0, MARGIN = 1, FUN = tangent_space_projection, B = G_avg))

  save(G_tangent, file=file.path(dir_features, 'G_tangent.RData'))

  # Run ComBat to harmonize U0 -------------------------------------------------

  # Project U0* back to get S0* ------------------------------------------------

  # Run ComBat to harmonize G_tangent ------------------------------------------

  # Transform G_tangent* back to get G0* ---------------------------------------

  # Rotate A0 so that Cov(A0*) = G0* -------------------------------------------

  # Reconstruct the fMRI data! -------------------------------------------------

  ### Compute the original residual E = Y - AS

  ### Compute empirically what % of variation in Y is left in E, return this too

  ### Compute Y* = A*S* + E

  # Organize stuff to return ---------------------------------------------------

  # harmonized Y
  # harmonized and unharmonized features?
  # distance between harmonized and unharmonized features?
  # var left in E
  # params (include GICA?)



  params <- list(
    #inds=inds,
    scale=scale,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    center_Bcols=center_Bcols,
    brainstructures=brainstructures,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol
  )

  list(DR=DR0, params=params)
}

library(expm)
tangent_space_projection <- function(A, B, reverse=FALSE) {
  # Assuming A and B are both positive definite matrices

  # Perform eigenvalue decomposition of B
  eig <- eigen(B)
  # Square root of eigenvalues
  sqrt_eigvals <- sqrt(eig$values)

  # Reconstruct square root of B
  sqrt_B <- eig$vectors %*% diag(sqrt_eigvals) %*% t(eig$vectors)
  inv_sqrt_B <- eig$vectors %*% diag(1/sqrt_eigvals) %*% t(eig$vectors)

  # Compute B^{-1/2}AB^{-1/2}
  middle_term <- inv_sqrt_B %*% A %*% inv_sqrt_B

  # Perform transformation or reverse transformation
  if(!reverse){
    # Compute the matrix logarithm
    log_term <- logm(middle_term)
    # Compute the projection
    projection <- sqrt_B %*% log_term %*% sqrt_B
  } else {
    exp_term <- expm(middle_term)
    projection <- sqrt_B %*% exp_term %*% sqrt_B
  }

  return(projection)
}

#' DR step for harmonize
#'
#' Compute dual regression for harmonize
#'
#' @param BOLD Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, GIFTI file paths,
#'  \code{"gifti"} objects, NIFTI file paths, \code{"nifti"} objects,
#'  or \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data
#'  locations and \eqn{T} is the number of timepoints.
#
#   If GIFTI or \code{"gifti"}, the input can also be a length two list,
#   where the first list element is a length \eqn{N} vector for the left
#   hemisphere and the second list element is a length \eqn{N} vector for the
#   right hemisphere.
#'  @param ghemi Which hemisphere, in the case of a single GIFTI file?
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one
#'  of the following: a \code{"CIFTI"} file path, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file path, a \code{"nifti"} object, or a \code{"data"} matrix.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also
#'  be a (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of
#'  \code{BOLD}. Its columns will be centered.
#' @param mask Required if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects, optional for other formats. For NIFTI, this is a brain map formatted as a
#'  binary array of the same spatial dimensions as the fMRI data, with
#'  \code{TRUE} corresponding to in-mask voxels. For other formats, a logical vector.
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation.
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents surface data (CIFTI or
#'  GIFTI). To smooth the standard deviation estimates used for local scaling,
#'  provide the surface geometries along which to smooth as GIFTI geometry files
#'  or \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
#'
#'  If \code{scale_sm_FWHM==0}, no smoothing of the local standard deviation
#'  estimates will be performed.
#'
#'  If \code{scale_sm_FWHM>0} but \code{scale_sm_surfL} and
#'  \code{scale_sm_surfR} are not provided, the default inflated surfaces from
#'  the HCP will be used.
#'
#'  To create a \code{"surf"} object from data, see
#'  \code{\link[ciftiTools]{make_surf}}. The surfaces must be in the same
#'  resolution as the \code{BOLD} data.
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
#' @param center_Bcols Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' @param maskTol For computing the dual regression results for each subject:
#'  tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  a subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for both scans and then
#'  the intersection of the masks is used, for each subject.
#' @param missingTol For computing the variance decomposition across all subjects:
#'  tolerance for number of subjects masked out due to low variance or missing
#'  values at a given location. If more than this many subjects are masked out,
#'  the location's value will be \code{NA} in the templates. \code{missingTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of subjects can be masked out
#'  at a given location.
#'
#' @keywords internal
harmonize_DR_oneBOLD <- function(
  BOLD, ghemi,
  format=c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"),
  GICA, mask=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  center_Bcols=FALSE,
  NA_limit=.1,
  brainstructures=c("all"),
  varTol=1e-6, maskTol=.1,
  verbose=TRUE){

  if (verbose) { extime <- Sys.time() }

  scale <- match.arg(scale, c("local", "global", "none"))
  # No other arg checks: check them before calling this function.

  # For `"xifti"` data for handling the medial wall and smoothing.
  xii1 <- NULL

  # Load helper variables.
  format <- match.arg(format, c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"))
  FORMAT <- fMRItools:::get_FORMAT(format)
  nQ <- ncol(GICA)
  nI <- nV <- nrow(GICA)

  fMRItools:::check_req_ifti_pkg(FORMAT)

  # Get `BOLD` (and `BOLD2`) as a data matrix or array.  -----------------------
  if (verbose) { cat("\tReading and formatting data...") }
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures) }
    if (ciftiTools::is.xifti(BOLD)) {
      if (scale == "local") {
        xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(BOLD, 1), "dscalar") * 0
      }
      BOLD <- as.matrix(BOLD)
    }
    stopifnot(is.matrix(BOLD))
    nI <- nV <- nrow(GICA)
  } else if (FORMAT == "GIFTI") {
    if (is.character(BOLD)) { BOLD <- gifti::readgii(BOLD) }
    stopifnot(gifti::is.gifti(BOLD))
    #ghemi <- BOLD$file_meta["AnatomicalStructurePrimary"]
    #if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
    #  stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    #}
    #ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    if (scale == "local") {
      if (ghemi == "left") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD$data)), 1) * 0
      } else if (ghemi == "right") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD$data)), 1) * 0
      } else { stop() }
      xii1$meta$cifti$intent <- 3006
    }
    BOLD <- do.call(cbind, BOLD$data)

    stopifnot(is.matrix(BOLD))
    nI <- nV <- nrow(GICA)
  } else if (FORMAT == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- RNifti::readNifti(BOLD) }
    stopifnot(length(dim(BOLD)) > 1)
    stopifnot(!is.null(mask))
    nI <- dim(drop(mask)); nV <- sum(mask)
  } else if (FORMAT == "MATRIX") {
    if (is.character(BOLD)) { BOLD <- readRDS(BOLD) }
    stopifnot(is.matrix(BOLD))
    nI <- nV <- nrow(GICA)
  } else { stop() }

  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))
  nT <- dim(BOLD)[ldB]

  # Check BOLD (and BOLD2) dimensions correspond with `GICA` and `mask`.
  if(!(ldB-1 == length(nI))) { stop("`GICA` and BOLD spatial dimensions do not match.") }
  if(!all(dBOLD[seq(ldB-1)] == nI)) { stop("`GICA` and BOLD spatial dimensions do not match.") }

  # Vectorize `BOLD`. ----------------------------------------------------------
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(as.logical(mask), dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
  }

  # Check for missing values. --------------------------------------------------
  nV0 <- nV # not used
  #mask <- fMRItools:::make_mask(BOLD, varTol=varTol)
  use_mask <- !all(mask)
  if (use_mask) {

    # # Coerce `maskTol` to number of locations.
    # stopifnot(is.numeric(maskTol) && length(maskTol)==1 && maskTol >= 0)
    # if (maskTol < 1) { maskTol <- maskTol * nV }
    # # Skip this scan if `maskTol` is surpassed.
    # if (sum(!mask) > maskTol) { return(NULL) }

    # Mask out the locations.
    BOLD <- BOLD[mask,,drop=FALSE]
    GICA <- GICA[mask,,drop=FALSE]
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask,] <- NA
      xii1 <- ciftiTools::move_to_mwall(ciftiTools::newdata_xifti(xii1, xiitmp))
    }
    nV <- nrow(BOLD)
  }

  if (!is.null(xii1) && scale=="local" && scale_sm_FWHM > 0) {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  # [TEMP? no scaling smoothing, no hpf] [TO DO]
  center_Bcols <- FALSE
  DR <- templateICAr::dual_reg(
    BOLD, GICA,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    #TR=TR, hpf=0,
    center_Bcols=center_Bcols
  )
  attr(DR$A, "scaled:center") <- NULL

  DR
}