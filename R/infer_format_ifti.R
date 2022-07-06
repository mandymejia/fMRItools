#' Infer fMRI data format
#'
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return The format: \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"NIFTI"} file path, \code{"nifti"} object, or \code{"data"}.
#' @keywords internal
infer_format_ifti <- function(BOLD, verbose=FALSE){

  # Character vector: CIFTI or NIFTI
  if (is.character(BOLD)) {
    stopifnot(length(BOLD)==1)
    if (endsWith(BOLD, ".dtseries.nii") | endsWith(BOLD, ".dscalar.nii")) {
      format <- "CIFTI"
    } else if (endsWith(BOLD, "gii")) {
      format <- "GIFTI"
    } else {
      format <- "NIFTI"
    }

  } else if (inherits(BOLD, "xifti")) {
    format <- "xifti"
  } else if (inherits(BOLD, "niftiExtension") | inherits(BOLD, "niftiImage") | inherits(BOLD, "nifti")) {
    format <- "nifti"
  } else if (inherits(BOLD, "gifti")) {
    format <- "gifti"

  } else {
    if (is.list(BOLD)) { stop("Unknown `BOLD` format.") }
    if (is.matrix(BOLD)) { 
      format <- "data"
    } else if (length(dim(BOLD))==4) {
      format <- "nifti"
    } else {
      stop("Unknown `BOLD` format.")
    }
  }
  if (verbose) { cat("Inferred input format:", format, "\n") }
  format
}