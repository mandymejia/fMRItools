#' Get NIFTI ROI masks
#'
#' Get NIFTI ROI masks
#'
#' @param nii_labels \eqn{I} by \eqn{J} by \eqn{K} NIFTI object or array (or
#'  file path to the NIFTI) which contains the corresponding labels to each
#'  voxel in \code{nii}. Values should be according to this table:
#'  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
#'  In the HCP, the corresponding file is "ROIs/Atlas_wmparc.2.nii.gz".
#' @param ROI_noise A list of numeric vectors. Each entry should represent
#'  labels in \code{nii_labels} belonging to a single noise ROI, named by that
#'  entry's name. Or, this can be a character vector of at least one of the
#'  following: \code{"wm_cort"} (cortical white matter), \code{"wm_cblm"}
#'  (cerebellar white matter), \code{"csf"} (cerebrospinal fluid). In the
#'  latter case, these labels will be used:
#'
#'  \describe{
#'    \item{\code{"wm_cort"}}{\code{c(3000:4035, 5001, 5002)}}
#'    \item{\code{"wm_cblm"}}{\code{c(7, 46)}}
#'    \item{\code{"csf"}}{\code{c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251,
#'      252, 253, 254, 255))}}
#'  }
#'
#'  These default ROIs are based on this forum post:
#'  https://www.mail-archive.com/hcp-users@humanconnectome.org/msg00931.html
#'
#'  Default: \code{c("wm_cort", "csf")}
#' @return The ROIs
#' @keywords internal
#'
get_NIFTI_ROI_masks <- function(nii_labels, ROI_noise=c("wm_cort", "csf")){
  # Read NIFTI labels file.
  if (is.character(nii_labels)) {
    nii_labels <- read_nifti(nii_labels)
  }
  stopifnot(length(dim(nii_labels))==3)

  # Get ROI masks.
  ROI_noise.default <- list(
    wm_cort = c(3000:4035, 5001, 5002),
    wm_cblm =  c(7, 46),
    csf = c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255)
  )
  if (is.null(ROI_noise)) {
    ROI_noise <- ROI_noise.default[c("wm_cort", "csf")]
  } else if (is.character(ROI_noise)) {
    stopifnot(all(ROI_noise %in% names(ROI_noise.default)))
    ROI_noise <- ROI_noise.default[unique(ROI_noise)]
  } else if (is.list(ROI_noise)) {
    ROI_noise <- ROI_noise
    stopifnot(is.numeric(do.call(c, ROI_noise)))
  } else {
    stop("`ROI_noise` argument is not in a recognized form.")
  }
  lapply(
    ROI_noise,
    function(x){array(nii_labels %in% x, dim=dim(nii_labels))}
  )
}

#' Anatomical CompCor for HCP NIFTI data
#'
#' Wrapper to \code{\link{CompCor}} for HCP-format data. Obtains the noise
#'  principal components (PCs) and ROIs from the NIFTI fMRI data and the NIFTI
#'  labels. To clean the surface-based CIFTI data with aCompCor using these
#'  noise PCs, see \code{\link{CompCor_HCP_cii}}.
#'
#' The NIFTI data and labels are read once, and the noise ROI masks (including
#'  erosion) and noise ROI data are extracted once. If \code{idx} is a list,
#'  the noise PCs are then computed separately for each entry of \code{idx}
#'  (e.g., for several scan durations) without re-reading or re-extracting the
#'  data.
#'
#' @inheritParams get_NIFTI_ROI_masks
#' @param nii \eqn{I} by \eqn{J} by \eqn{K} by \eqn{T} NIFTI object or array
#'  (or file path to the NIFTI) which contains whole-brain data, including the
#'  noise ROIs. In the HCP, the corresponding file is e.g.
#'  "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz"
#' @param noise_nPC The number of principal components to compute for each
#'  noise ROI. Alternatively, values between 0 and 1, in which case they will
#'  represent the minimum proportion of variance explained by the PCs used for
#'  each noise ROI. The smallest number of PCs will be used to achieve this
#'  proportion of variance explained.
#'
#'  Should be a list or numeric vector with the same length as
#'  \code{ROI_noise}. It will be matched to each ROI based on the name of each
#'  entry, or if the names are missing, the order of entries. If it is an
#'  unnamed vector, its elements will be recycled. Default: \code{5} (compute
#'  the top 5 PCs for each noise ROI).
#' @param noise_erosion  The number of voxel layers to erode the noise ROIs by.
#'  Should be a list or numeric vector with the same length as
#'  \code{ROI_noise}. It will be matched to each ROI based on the name of each
#'  entry, or if the names are missing, the order of entries. If it is an
#'  unnamed vector, its elements will be recycled. Default: \code{NULL}, which
#'  will use a value of 0 (do not erode the noise ROIs).
#' @param idx A numeric vector indicating the timepoints to use, or
#'  \code{NULL} (default) to use all idx. (Indexing begins with 1, so the
#'  first timepoint has index 1 and the last has the same index as the length
#'  of the scan.)
#'
#'  Alternatively, a list of such numeric vectors (e.g. one for each scan
#'  duration of interest). In this case the noise PCs are computed separately
#'  for each entry (reading the data only once) and a list with one result per
#'  entry is returned, with the same names as \code{idx}.
#' @param center,scale Center the columns of the noise ROI data by median, and
#'  scale the columns of the noise ROI data by MAD? Default: \code{TRUE} for
#'  both.
#' @param prelim_nuisance Nuisance signals to regress from the noise ROIs
#'  before computing the PCs (see \code{\link{CompCor}}). A \eqn{T} by
#'  \eqn{P} numeric matrix, where \eqn{T} is the number of timepoints in
#'  \code{idx}. If \code{idx} is a list, this must be a list of matrices of
#'  the same length, with the \eqn{k}th matrix having as many rows as the
#'  \eqn{k}th entry of \code{idx}. Default: \code{NULL} (do not do any
#'  preliminary nuisance regression).
#' @param verbose Should occasional updates be printed? Default:
#'  \code{FALSE}.
#'
#' @return A list with entries \code{"PCs"} (the noise PC scores for each
#'  noise ROI), \code{"var"} (the variance of each PC), and \code{"ROI_noise"}
#'  (the noise ROI masks, after erosion). If \code{idx} is a list, a list of
#'  these, one for each entry of \code{idx}.
#'
#' @section References:
#'  \itemize{
#'    \item{Behzadi, Y., Restom, K., Liau, J. & Liu, T. T. A component based
#'      noise correction method (CompCor) for BOLD and perfusion based fMRI.
#'      NeuroImage 37, 90-101 (2007).}
#'    \item{Muschelli, J. et al. Reduction of motion-related artifacts in
#'      resting state fMRI using aCompCor. NeuroImage 96, 22-35 (2014).}
#' }
#'
#' @export
#'
#' @seealso CompCor, CompCor_HCP_cii
CompCor_HCP <- function(
  nii, nii_labels,
  ROI_noise=c("wm_cort", "csf"), noise_nPC=5, noise_erosion=NULL,
  idx=NULL, center = TRUE, scale = TRUE,
  prelim_nuisance = NULL, verbose=FALSE){

  # Get NIFTI data (read only once).
  if (is.character(nii)) {
    if (verbose) { cat("Reading data NIFTI.\n") }
    nii <- read_nifti(nii)
  }
  stopifnot(length(dim(nii))==4)
  T_nii <- dim(nii)[4]

  # Get NIFTI ROI masks.
  if (verbose) { cat("Reading labels NIFTI.\n") }
  ROI_noise <- get_NIFTI_ROI_masks(nii_labels, ROI_noise)
  stopifnot(all( dim(ROI_noise[[1]]) == dim(nii)[1:3] ))

  # Standardize `idx` to a list of index vectors.
  multi_idx <- is.list(idx)
  if (is.null(idx)) {
    idx <- list(seq(T_nii))
  } else if (!multi_idx) {
    idx <- list(idx)
  }
  for (kk in seq_along(idx)) {
    stopifnot(is.numeric(idx[[kk]]))
    stopifnot(all(idx[[kk]] %in% seq(T_nii)))
    if (length(idx[[kk]]) < 10) {
      warning("There are very few `idx`.\n")
    }
  }

  # Standardize `prelim_nuisance` to a list with one entry per `idx` entry.
  if (is.null(prelim_nuisance)) {
    prelim_nuisance <- vector("list", length(idx))
  } else {
    if (!is.list(prelim_nuisance) || is.data.frame(prelim_nuisance)) {
      prelim_nuisance <- list(as.matrix(prelim_nuisance))
    }
    if (length(prelim_nuisance) != length(idx)) {
      stop(
        "`prelim_nuisance` must be a matrix, or a list with one matrix per ",
        "entry of `idx`."
      )
    }
    for (kk in seq_along(idx)) {
      if (nrow(prelim_nuisance[[kk]]) != length(idx[[kk]])) {
        stop(
          "`prelim_nuisance` entry ", kk, " does not have the same number of ",
          "rows as the length of `idx` entry ", kk, "."
        )
      }
    }
  }

  # Get the noise ROI data (with erosion) once, for all timepoints.
  if (verbose) { cat("Extracting noise ROIs.\n") }
  fmt <- format_data(
    X=nii, ROI_data=NULL, ROI_noise=ROI_noise,
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )
  rm(nii)

  # Compute the noise components for each `idx`.
  noise <- vector("list", length(idx))
  names(noise) <- names(idx)
  for (kk in seq_along(idx)) {
    if (verbose) {
      cat("Computing noise components (", kk, " of ", length(idx), ").\n",
        sep="")
    }

    # Subset timepoints.
    X_noise <- lapply(fmt$X_noise, function(x){ x[idx[[kk]], , drop=FALSE] })

    # Preliminary nuisance regression. This is the same projection as
    #   `nuisance_regression` (rank-deficient designs are handled in the same
    #   way as in `hat_matrix`), but it uses the thin QR factor of the design
    #   matrix instead of the T x T hat matrix, and only computes it once for
    #   all ROIs.
    if (!is.null(prelim_nuisance[[kk]])) {
      design <- as.matrix(prelim_nuisance[[kk]])
      qrd <- qr(design)
      design <- design[, qrd$pivot[seq_len(qrd$rank)], drop=FALSE]
      Qd <- qr.Q(qr(design))
      X_noise <- lapply(X_noise, function(x){ x - Qd %*% crossprod(Qd, x) })
    }

    out <- CompCor.noise_comps(
      X_noise=X_noise, center=center, scale=scale, noise_nPC=fmt$noise_nPC
    )
    noise[[kk]] <- list(
      PCs=out$noise_comps, var=out$noise_var, ROI_noise=fmt$ROI_noise
    )
  }

  if (multi_idx) { noise } else { noise[[1]] }
}

#' Clean HCP CIFTI data with aCompCor
#'
#' Regresses noise PCs from HCP CIFTI data. The PCs can be calculated from the
#'  NIFTI fMRI data and NIFTI labels with \code{\link{CompCor_HCP}}.
#'
#' @param cii \code{"xifti"} (or file path to the CIFTI) from which the noise
#'  ROI components will be regressed. In the HCP, the corresponding file is
#'  e.g.
#'  "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii".
#' @param noise The noise components, as returned by \code{\link{CompCor_HCP}}
#'  for a single \code{idx} (i.e. a list with an entry named \code{"PCs"}).
#' @param idx A numeric vector indicating the timepoints of \code{cii} to use,
#'  or \code{NULL} (default) to use all timepoints. This should match the
#'  timepoints used to compute \code{noise}.
#' @param brainstructures Choose among "left", "right", and "subcortical".
#'  Only used if \code{cii} is a file path. Default: \code{c("left", "right")}
#'  (cortical data only)
#' @param DCT Add DCT bases to the nuisance regression? Use an integer to
#'  indicate the number of cosine bases. Use \code{0} (default) to forgo
#'  detrending.
#'
#'  The data must be centered, either before input or with \code{center}.
#' @param nuisance_too A matrix of nuisance signals to add to the nuisance
#'  regression. Should have \eqn{T} rows. \code{NULL} to not add additional
#'  nuisance regressors (default).
#' @param verbose Should occasional updates be printed? Default:
#'  \code{FALSE}.
#'
#' @return The cleaned surface-based data as a \code{"xifti"} object.
#'
#' @export
#'
#' @seealso CompCor_HCP
CompCor_HCP_cii <- function(
  cii, noise, idx=NULL, brainstructures=c("left", "right"),
  DCT = 0, nuisance_too = NULL, verbose=FALSE){

  stopifnot(is.list(noise), !is.null(noise$PCs))
  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package `ciftiTools` is required. Please install it.")
  }

  # Get CIFTI data.
  if (is.character(cii)) {
    cii <- ciftiTools::read_cifti(
      cii, brainstructures=brainstructures, verbose=verbose
    )
  }
  stopifnot(all(names(cii) == c("data", "surf", "meta")))

  # Drop idx.
  if (!is.null(idx)) { cii <- ciftiTools::select_xifti(cii, idx) }
  cii_mat <- as.matrix(cii)
  T_ <- ncol(cii_mat)

  # Make design matrix.
  design <- do.call(cbind, noise$PCs)
  if (!is.null(design) && nrow(design) != T_) {
    stop(
      "The noise PCs have ", nrow(design), " timepoints, but `cii` has ", T_,
      " (after applying `idx`)."
    )
  }
  if (DCT > 0) { design <- cbind(design, dct_bases(T_, DCT)) }
  if (!is.null(nuisance_too)) {
    nuisance_too <- validate_design_mat(nuisance_too, T_)
    design <- cbind(design, nuisance_too)
  }
  design <- validate_design_mat(cbind(1, design), T_)

  # Return data in "xifti" format.
  ciftiTools::newdata_xifti(cii, nuisance_regression(cii_mat, design))
}
