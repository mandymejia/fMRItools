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
#'  locations and \eqn{T} is the number of timepoints. For GIFTI or
#'  \code{"gifti"} input, each entry can also be a length-two list where the
#'  first entry is the left cortex and the second is the right cortex.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also
#'  be a (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of
#'  \code{BOLD}. Its columns will be centered.
#' @param mask Required if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects, optional for other formats. For
#'  NIFTI, this is a brain map formatted as a binary array of the same spatial
#'  dimensions as the fMRI data, with \code{TRUE} corresponding to in-mask
#'  voxels. For other formats, a logical vector.
#' @param gii_hemi Which hemisphere, if BOLD is a single list of GIFTI data?
#'  Should be \code{"left"} or \code{"right"}. If \code{NULL} (default), try to
#'  infer from the metadata.
#' @param inds Numeric indices of the group ICs to include in the template. If
#'  \code{NULL}, use all group ICs (default).
#'
#'  If \code{inds} is provided, the ICs not included will be removed after
#'  calculating dual regression, not before. This is because removing the ICs
#'  prior to dual regression would leave unmodelled signals in the data, which
#'  could bias the templates.
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
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' @param do_harmonize Perform the harmonization?  If not, will only extract and
#' return features to be harmonized.
#' @keywords internal
#' @importFrom stats cov
#'
harmonize <- function(
  BOLD, GICA,
  mask=NULL,
  gii_hemi=NULL,
  inds=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL,
  scale_sm_surfR=NULL,
  scale_sm_FWHM=2,
  TR=NULL,
  hpf=.01,
  GSR=FALSE,
  brainstructures=c("all"),
  verbose=TRUE,
  do_harmonize=FALSE){

  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("Package \"expm\" needed. Please install it", call. = FALSE)
  }
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
  stopifnot(is_1(scale_sm_FWHM, "numeric"))
  if (is.null(hpf)) { hpf <- 0 }
  if (is.null(TR)) {
    if (hpf==.01) {
      message("Setting `hpf=0` because `TR` was not provided. Either provide `TR` or set `hpf=0` to disable this message.")
      hpf <- 0
    } else {
      stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=0`.")
    }
  } else {
    stopifnot(is_posNum(TR))
    stopifnot(is_posNum(hpf, zero_ok=TRUE))
  }
  stopifnot(is_1(GSR, "logical"))
  stopifnot(is_1(verbose, "logical"))

  # `BOLD` format --------------------------------------------------------------
  format <- infer_format_ifti_vec(BOLD)[1]
  FORMAT <- get_FORMAT(format)
  FORMAT_extn <- switch(FORMAT,
    CIFTI=".dscalar.nii",
    GIFTI=".func.gii",
    GIFTI2=".func.gii",
    NIFTI=".nii",
    MATRIX=".rds"
  )
  nN <- length(BOLD)

  check_req_ifti_pkg(FORMAT)

  # If BOLD is a CIFTI, GIFTI, GIFTI2, NIFTI, or RDS file, check that the file paths exist.
  if (format %in% c("CIFTI", "GIFTI", "GIFTI2", "NIFTI", "RDS")) {

    if (format == "GIFTI2") {
      missing_BOLD1 <- !file.exists(vapply(BOLD, function(q){ q[[1]] }, ''))
      missing_BOLD2 <- !file.exists(vapply(BOLD, function(q){ q[[2]] }, ''))
      missing_BOLD <- missing_BOLD1 | missing_BOLD2
    } else {
      missing_BOLD <- !file.exists(BOLD)
    }

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
    gii_hemi <- GICA$file_meta["AnatomicalStructurePrimary"]
    if (!(gii_hemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    }
    gii_hemi <- switch(gii_hemi, CortexLeft="left", CortexRight="right")
    GICA <- do.call(cbind, GICA$data)
  } else if (FORMAT == "GIFTI2") {
    GICA <- as.list(GICA)
    if (length(GICA) != 2) {
      stop("`GICA` should be a length-2 list of GIFTI data: left hemisphere first, right hemisphere second.")
    }
    for (ii in seq(2)) {
      if (is.character(GICA[[ii]])) { GICA[[ii]] <- gifti::readgii(GICA[[ii]]) }
      gii_hemi_ii <- GICA[[ii]]$file_meta["AnatomicalStructurePrimary"]
      if (!(gii_hemi_ii %in% c("CortexLeft", "CortexRight"))) {
        stop("AnatomicalStructurePrimary metadata missing or invalid for GICA[[ii]].")
      }
      gii_hemi_ii <- switch(gii_hemi_ii, CortexLeft="left", CortexRight="right")
      if (gii_hemi_ii != c("left", "right")[ii]) {
        stop("`GICA` should be a length-2 list of GIFTI data: left hemisphere first, right hemisphere second.")
      }
      GICA[[ii]] <- do.call(cbind, GICA[[ii]]$data)
    }
    GICA <- rbind(GICA[[1]], GICA[[2]])
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- RNifti::readNifti(GICA) }
    stopifnot(length(dim(GICA)) > 1)
  } else if (FORMAT == "MATRIX") {
    if (is.character(GICA)) { GICA <- readRDS(GICA) }
    stopifnot(is.matrix(GICA))
  }
  nQ <- dim(GICA)[length(dim(GICA))]

  # `inds`.
  if (!is.null(inds)) {
    if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
    nL <- length(inds)
  } else {
    inds <- seq(nQ)
    nL <- nQ
  }

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
  } else { #For non-NIFTI data, mask is not required but can be provided
    if (!is.null(mask)) {
      #OPTIONAL GIFTI-FORMAT MASK
      if (FORMAT == "GIFTI") {
        # DAMON FILL IN HERE
      } else if (FORMAT == "GIFTI2") {
        mask <- as.list(mask)
        if (length(mask) != 2) {
          stop("`mask` should be a length-2 list of GIFTI data: left hemisphere first, right hemisphere second.")
        }
        for (ii in seq(2)) {
          if (is.character(mask[[ii]])) { mask[[ii]] <- gifti::readgii(mask[[ii]]) }
          gii_hemi_ii <- mask[[ii]]$file_meta["AnatomicalStructurePrimary"]
          if (!(gii_hemi_ii %in% c("CortexLeft", "CortexRight"))) {
            stop("AnatomicalStructurePrimary metadata missing or invalid for mask[[ii]].")
          }
          gii_hemi_ii <- switch(gii_hemi_ii, CortexLeft="left", CortexRight="right")
          if (gii_hemi_ii != c("left", "right")[ii]) {
            stop("`mask` should be a length-2 list of GIFTI data: left hemisphere first, right hemisphere second.")
          }
          mask[[ii]] <- do.call(cbind, mask[[ii]]$data)
        }
        mask <- rbind(mask[[1]], mask[[2]])
        mask <- as.logical(mask)
      } else if (FORMAT == "MATRIX") {
        # DAMON FILL IN HERE
      }
      nI <- length(mask); nV <- sum(mask)
    } else { #end if(!is.null(mask))
      nI <- nV <- nrow(GICA)
    }
  } #end else (not NIFTI format)

  # [TO DO]: Apply mask to GICA prior to centering

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- colCenter(GICA) }

  # Print summary of data ------------------------------------------------------
  format2 <- if (format == "data") { "numeric matrix" } else { format }
  if (verbose) {
    cat("Data input format:    ", format2, "\n")
    cat("Image dimensions:     ", paste(nI, collapse=" x "), "\n")
    cat("In-mask locations:    ", nV, "\n")
    if (FORMAT == "GIFTI") {
      cat("Cortex hemisphere:    ", gii_hemi, "\n")
    }
    cat("Number of group ICs:  ", nQ, "\n")
    cat("Number of sessions:   ", nN, "\n")
  }

  # [TO DO]: Implement template ICA instead of DR. The estimates are so noisy.

  # Dual regression ------------------------------------------------------------
  S0 <- array(0, dim = c(nN, nQ, nV))
  A0 <- vector("list", nN)
  G0 <- array(NA, dim = c(nN, nQ, nQ))
  for (ii in seq(nN)) {
    if (verbose) { cat(paste0(
      '\nReading and analyzing data for subject ', ii,' of ', nN, '.\n'
    )) }

    DR0_ii <- harmonize_DR_oneBOLD(
      BOLD[[ii]],
      mask=mask,
      gii_hemi=gii_hemi,
      format=format,
      GICA=GICA,
      GSR=GSR,
      scale=scale,
      scale_sm_surfL=scale_sm_surfL,
      scale_sm_surfR=scale_sm_surfR,
      scale_sm_FWHM=scale_sm_FWHM,
      TR=TR,
      hpf=hpf,
      brainstructures=brainstructures,
      verbose=verbose
    )
    S0[ii,,] <- DR0_ii$S[,mask]
    A0[[ii]] <- DR0_ii$A
    G0[ii,,] <- cov(DR0_ii$A)
  }

  ### Process the features -------------------------------------------------------

    ## Use SVD to reduce dimensions of S_q <- MAYBE NOT!  USE SPATIAL MODEL INSTEAD?

  # U0 <- V0 <- vector('list', nQ)
  # for (qq in seq(nQ)) {
  #
  #   if (verbose) { cat(paste0(
  #     '\nPerforming PCA on IC ', qq,' of ', nQ, '.\n'
  #   )) }
  #
  #   #do PCA separately for each IC q
  #   S0_q <- S0[,qq,] #NxV
  #
  #   #center across sessions so X'X is covariance matrix (add back later)
  #   #S0_q_mean <- colMeans(S0_q)
  #   S0_q <- t(S0_q) #- S0_q_mean #vector will be recycled by column
  #
  #   #want to obtain U (NxP) where P << V
  #   SSt_q <- crossprod(S0_q) # S0 is currently VxN, so S0'S0 is NxN.
  #   svd_q <- svd(SSt_q, nv=0) #want to get U from SVD of S0' = UDV'. SVD of S0'S0 = UDV'VDU' = U D^2 U'
  #
  #   #keep components explaining 99% of variation
  #   nP <- min(which(cumsum(svd_q$d)/sum(svd_q$d) > 0.99))
  #   U0[[qq]] <- svd_q$u[,1:nP]
  #
  #   #visualize V' = (1/D)U'S0' to make sure they are sensible
  #   V0[[qq]] <- diag(1/svd_q$d[1:nP]) %*% t(svd_q$u[,1:nP]) %*% t(S0_q)
  # }

  ## Project elements of G to a tangent space

  if (verbose) { cat('\nProjecting covariance matrices to tangent space.\n') }

  # Calculate the element-wise average of the covariance matrices
  G_avg <- apply(G0, c(2, 3), mean)

  # Project each covariance matrix to tangent space at base point defined by the Euclidean mean
  G_tangent <- apply(G0, MARGIN = 1, FUN = tangent_space_projection, B = G_avg) #very fast

  feature_list <- list(S = S0,
                       A = A0,
                       G = G0,
                       Gt = G_tangent)

  if(!do_harmonize) return(feature_list)


  # Run ComBat to harmonize S0 -------------------------------------------------

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
    inds=inds,
    scale=scale,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    GSR=GSR,
    brainstructures=brainstructures
  )

  # Harmonize A and cov(A) using ComBat — given the S and cov(A) matrices, Johanna can develop the code to do that
  # Re-construct the fMRI data as Y* = A*S* + E, where A* and S* are the harmonized versions of A and S.  The harmonized version of A is done using matrix operations so that cov(A*) = cov(A)* (the harmonized covariant matrix).  — you could write this step as well, and just leave a placeholder for the harmonization itself
  # Return/write out the harmonized fMRI time series

  list(
    features=feature_list,
    params=params
  )
}

#' Tangent space projection
#'
#' Tangent space projection
#'
#' @param A,B,reverse To-Do
#' @return To-Do
#' @keywords internal
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
    log_term <- expm::logm(middle_term)
    # Compute the projection
    projection <- sqrt_B %*% log_term %*% sqrt_B
  } else {
    exp_term <- expm::expm(middle_term)
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
#'  @param gii_hemi Which hemisphere, in the case of a single GIFTI file?
#' @param format Expected format of \code{BOLD}. Should be one
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
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#'
#' @keywords internal
harmonize_DR_oneBOLD <- function(
  BOLD,
  format=c("CIFTI", "xifti", "GIFTI", "gifti", "GIFTI2", "gifti2", "NIFTI", "nifti", "RDS", "data"),
  GICA,
  mask=NULL,
  gii_hemi=NULL,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE,
  NA_limit=.1,
  brainstructures=c("all"),
  verbose=TRUE){

  if (verbose) { extime <- Sys.time() }

  if (!is.null(gii_hemi)) { stopifnot(gii_hemi %in% c("left", "right")) }

  scale <- match.arg(scale, c("local", "global", "none"))
  # No other arg checks: check them before calling this function.

  # For handling cortical surface data: the medial wall and smoothing.
  xii1 <- NULL

  # Load helper variables.
  format <- match.arg(format, c("CIFTI", "xifti", "GIFTI", "gifti", "GIFTI2", "gifti2", "NIFTI", "nifti", "RDS", "data"))
  FORMAT <- get_FORMAT(format)
  check_req_ifti_pkg(FORMAT)

  nQ <- ncol(GICA)

  if(is.null(mask)){
    nI <- nV <- nrow(GICA)
  } else {
    nI <- nrow(GICA); nV <- sum(mask)
  }

  # Get `BOLD` as a data matrix or array.  -------------------------------------
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
    #nI <- nV <- nrow(GICA)
  } else if (FORMAT == "GIFTI") {
    if (is.character(BOLD)) { BOLD <- gifti::readgii(BOLD) }
    stopifnot(gifti::is.gifti(BOLD))
    if (is.null(gii_hemi)) {
      gii_hemi <- BOLD$file_meta["AnatomicalStructurePrimary"]
      if (!(gii_hemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
      }
      gii_hemi <- switch(gii_hemi, CortexLeft="left", CortexRight="right")
    } else {
      if (gii_hemi == "left" && BOLD$file_meta["AnatomicalStructurePrimary"]=="CortexRight") {
        stop("`gii_hemi` is 'left' but GIFTI data is for the right cortex.")
      }
      if (gii_hemi == "right" && BOLD$file_meta["AnatomicalStructurePrimary"]=="CortexLeft") {
        stop("`gii_hemi` is 'right' but GIFTI data is for the left cortex.")
      }
    }
    if (scale == "local") {
      if (gii_hemi == "left") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD$data)), 1) * 0
      } else if (gii_hemi == "right") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD$data)), 1) * 0
      } else { stop() }
      xii1$meta$cifti$intent <- 3006
    }
    BOLD <- do.call(cbind, BOLD$data)
    stopifnot(is.matrix(BOLD))
    #nI <- nV <- nrow(GICA)
  } else if (FORMAT == "GIFTI2") {
    BOLD_new <- vector('list', length=2)
    for (ii in seq(2)) {
      if (is.character(BOLD[[ii]])) { BOLD_new[[ii]] <- gifti::readgii(BOLD[[ii]]) }
      stopifnot(gifti::is.gifti(BOLD_new[[ii]]))
    }
    BOLD <- lapply(BOLD_new, function(q){do.call(cbind, q$data)})
    if (scale == "local") {
      xii1 <- ciftiTools::select_xifti(
        ciftiTools::as.xifti(
          cortexL=BOLD[[1]],
          cortexR=BOLD[[2]]
        ), 1
      ) * 0
      xii1$meta$cifti$intent <- 3006
      # the medial wall is included in the data.
      xii1$meta$cortex$resamp_res <- nrow(xii1$data$cortex_left)
    }
    BOLD <- do.call(rbind, BOLD)
  } else if (FORMAT == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- RNifti::readNifti(BOLD) }
    stopifnot(length(dim(BOLD)) > 1)
    stopifnot(!is.null(mask))
    nI <- dim(drop(mask)) #; nV <- sum(mask)
  } else if (FORMAT == "MATRIX") {
    if (is.character(BOLD)) { BOLD <- readRDS(BOLD) }
    stopifnot(is.matrix(BOLD))
    #nI <- nV <- nrow(GICA)
  } else { stop() }

  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))
  nT <- dim(BOLD)[ldB]

  # Check BOLD dimensions correspond with `GICA` and `mask`.
  if(!(ldB-1 == length(nI))) { stop("`GICA` and BOLD spatial dimensions do not match.") }
  if(!all(dBOLD[seq(ldB-1)] == nI)) { stop("`GICA` and BOLD spatial dimensions do not match.") }

  # Vectorize `BOLD`. ----------------------------------------------------------
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(as.logical(mask), dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
  }

  # Check for missing values. --------------------------------------------------
  mask <- as.logical(mask)
  use_mask <- !all(mask)
  if (use_mask) {
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

  GSR <- FALSE
  DR <- dual_reg(
    BOLD, GICA,
    scale=scale,
    scale_sm_xifti=xii1,
    scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf,
    GSR=GSR
  )
  attr(DR$A, "scaled:center") <- NULL

  if (use_mask) { DR$S <- unmask_mat(DR$S, mask, mask_dim=2, fill=NA) }

  DR
}
