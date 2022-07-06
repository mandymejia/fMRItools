#' Get NIFTI ROI masks
#'
#' Get NIFTI ROI masks
#'
#' @param nii_labels \eqn{I} by \eqn{J} by \eqn{K}
#'  NIFTI object or array (or file path to the NIFTI) which
#'  contains the corresponding labels to each voxel in \code{nii}. Values should
#'  be according to this table: 
#'  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT .
#'  In the HCP, the corresponding file is "ROIs/Atlas_wmparc.2.nii.gz". 
#' @param ROI_noise A list of numeric vectors. Each entry should represent labels
#'  in \code{nii_labels} belonging to a single noise ROI, named by that entry's 
#'  name. Or, this can be a character vector of at least one of the following: 
#'  \code{"wm_cort"} (cortical white matter), \code{"wm_cblm"} (cerebellar white
#'  matter), \code{"csf"} (cerebrospinal fluid). In the latter case, these labels
#'  will be used:
#'
#'  \describe{
#'    \item{\code{"wm_cort"}}{\code{c(3000:4035, 5001, 5002)}}
#'    \item{\code{"wm_cblm"}}{\code{c(7, 46)}}
#'    \item{\code{"csf"}}{\code{c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 251, 252, 253, 254, 255))}}
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

#' Anatomical CompCor for HCP NIFTI and CIFTI data
#'
#' Wrapper to \code{\link{CompCor}} for HCP-format data. Can be used to clean
#'  the surface-based CIFTI data with aCompCor using the noise PCs and ROIs 
#'  calculated from the NIFTI fMRI data and NIFTI mask. Can also be used to just
#'  obtain the noise PCs and ROIs without performing aCompCor, if the CIFTI
#'  data is not provided.
#' 
#' @inheritParams get_NIFTI_ROI_masks
#' @param nii \eqn{I} by \eqn{J} by \eqn{K} by \eqn{T} 
#'  NIFTI object or array (or file path to the NIFTI) which contains
#'  whole-brain data, including the noise ROIs. In the HCP, the corresponding
#'  file is e.g. "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz"
#' @param noise_nPC The number of principal components to compute for each noise
#'  ROI. Alternatively, values between 0 and 1, in which case they will 
#'  represent the minimum proportion of variance explained by the PCs used for
#'  each noise ROI. The smallest number of PCs will be used to achieve this 
#'  proportion of variance explained. 
#' 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its
#'  elements will be recycled. Default: \code{5} (compute the top 5 PCs for 
#'  each noise ROI).
#' @param noise_erosion  The number of voxel layers to erode the noise ROIs by. 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its 
#'  elements will be recycled. Default: \code{NULL}, which will use a value of
#'  0 (do not erode the noise ROIs).
#' @param brainstructures Choose among "left", "right", and "subcortical".
#'  Default: \code{c("left", "right")} (cortical data only)
#' @param idx A numeric vector indicating the timepoints to use, or 
#'  \code{NULL} (default) to use all idx. (Indexing begins with 1, so the 
#'  first timepoint has index 1 and the last has the same index as the length of 
#'  the scan.)
#' @param cii \code{"xifti"} (or file path to the CIFTI) from which the noise
#'  ROI components will be regressed. In the HCP, the corresponding file is e.g.
#'  "../Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll.dtseries.nii". If not
#'  provided, only the noise components will be returned (no data will be cleaned).
#' @param center,scale Center the columns of the data by median, and scale the
#'  columns of the data by MAD? Default: \code{TRUE} for both. Affects both
#'  \code{X} and the noise data. \code{center} also applies to \code{nuisance_too}
#'  so if it is \code{FALSE}, \code{nuisance_too} must already be centered.
#' @param DCT Add DCT bases to the nuisance regression? Use an integer to 
#'  indicate the number of cosine bases. Use \code{0} (default) to forgo detrending. 
#' 
#'  The data must be centered, either before input or with \code{center}.
#' @param nuisance_too A matrix of nuisance signals to add to the nuisance
#'  regression. Should have \eqn{T} rows. \code{NULL} to not add additional 
#'  nuisance regressors (default).
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @return The noise components, and if \code{cii} is provided, the cleaned
#'  surface-based data as a \code{"xifti"} object.
#' 
#' @section References:
#'  \itemize{
#'    \item{Behzadi, Y., Restom, K., Liau, J. & Liu, T. T. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. NeuroImage 37, 90-101 (2007).}
#'    \item{Muschelli, J. et al. Reduction of motion-related artifacts in resting state fMRI using aCompCor. NeuroImage 96, 22-35 (2014).}
#' }
#' 
#' @export
#' 
#' @seealso CompCor
CompCor_HCP <- function(
  nii, nii_labels, 
  ROI_noise=c("wm_cort", "csf"), noise_nPC=5, noise_erosion=NULL, 
  idx=NULL, cii=NULL, brainstructures=c("left", "right"),
  center = TRUE, scale = TRUE, DCT = 0, nuisance_too = NULL,
  verbose=FALSE){

  # Get NIFTI data.
  if (is.character(nii)) {
    if (verbose) { cat("Reading data NIFTI.\n") }
    nii <- read_nifti(nii)
  }
  stopifnot(length(dim(nii))==4)

  # Get NIFTI ROI masks.
  if (verbose) { cat("Reading labels NIFTI.\n") }
  ROI_noise <- get_NIFTI_ROI_masks(nii_labels, ROI_noise)
  stopifnot(all( dim(ROI_noise[[1]]) == dim(nii)[1:3] ))

  # Drop idx in NIFTI.
  T_ <- dim(nii)[4]
  if (!is.null(idx)) {
    stopifnot(all(idx %in% seq(T_)))
    nii <- nii[,,,idx]
  }
  T_ <- dim(nii)[4]
  if (T_ < 10) {
    warning("There are very few `idx`.\n")
  }

  if (verbose) { cat("Computing noise components.\n") }
  out <- CompCor(
    nii, ROI_data=NULL, ROI_noise=ROI_noise, 
    noise_erosion=noise_erosion, noise_nPC=noise_nPC,
    center=center, scale=scale,
    # We do nuisance regression with aCompCor, not prior.
    nuisance=NULL
  )

  if (is.null(cii)) { return(out$noise) }

  if (is.character(cii)) { 
    if (requireNamespace("ciftiTools", quietly = TRUE)) {
      cii <- ciftiTools::read_cifti(cii, brainstructures=brainstructures, verbose=verbose) 
    } else {
      stop("Package `ciftiTools` required to read the CIFTI file. Please install it.")
    }
  }
  stopifnot(all(names(cii) == c("data", "surf", "meta")))

  # Drop idx.
  if (!is.null(idx)) { cii <- ciftiTools::select_xifti(cii, idx) }

  # Make design matrix.
  design <- do.call(cbind, out$noise$PCs)
  if (DCT > 0) { design <- cbind(design, dct_bases(T_, DCT)) }
  if (!is.null(nuisance_too)) {
    nuisance_too <- check_design_matrix(nuisance_too)
    design <- cbind(design, nuisance_too)
  }
  design <- check_design_matrix(cbind(1, design))

  # Return data in \code{"xifti"} format.
  out$data <- ciftiTools::newdata_xifti(cii, nuisance_regression(cii, design))
  out
}