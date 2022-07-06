#' Wrapper to common functions for reading NIFTIs
#' 
#' Tries \code{RNifti::readNifti}, then \code{oro.nifti::readNIfTI}. If
#'  neither package is available an error is raised.
#' 
#' @param nifti_fname The file name of the NIFTI.
#' @return The NIFTI
#' @keywords internal
read_nifti <- function(nifti_fname){
  if (requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::readNifti(nifti_fname))
  } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
    return(oro.nifti::readNIfTI(nifti_fname, reorient=FALSE))
  } else {
    stop(
      "Package \"RNifti\" or \"oro.nifti\" needed to read", nifti_fname, 
      ". Please install at least one", call. = FALSE
    )
  }
}