#' Infer fMRI data format
#'
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return The format: \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"GIFTI"} file path, \code{"gifti"} object,  
#'  \code{"NIFTI"} file path, \code{"nifti"} object, 
#'  or \code{"data"}.
#' @keywords internal
infer_format_ifti <- function(BOLD, verbose=FALSE){

  Bformat <- Bformat2 <- NA

  # Define some file extensions to look out for
  cii_ext_supported <- paste0(".", ciftiTools::supported_intents()$extension)
  cii_ext_notsupport <- paste0(
    ".", c(
      "dconn", "pconn", "ptseries", "dtraj", "pscalar", 
      "pdconn", "dpconn", "pconnseries", "pconnscalar"
    ), ".nii"
  )
  gii_ext_supported <- c("func", "label")
  nii_classes <- c(oro.nifti="niftiExtension", RNifti="niftiImage", oro.nifti="nifti")

  # Character vector: CIFTI, GIFTI, or NIFTI
  if (is.character(BOLD) && length(BOLD)==1) {
    if (any(endsWith(BOLD, cii_ext_supported))) {
      for (cc in cii_ext_supported) {
        if (endsWith(BOLD, cc)) {
          Bformat <- "CIFTI"
          Bformat2 <- gsub("\\.nii$", "", gsub("^\\.", "", cc))
        }
      }
    } else if (any(endsWith(BOLD, cii_ext_notsupport))) {
      stop(
        "`BOLD` is a CIFTI file but only dtseries, dscalar, and dlabel ",
        "intents are supported."
      )
    } else if (endsWith(BOLD, ".gii")) {
      Bformat <- "GIFTI"
      gii_ext_supported <- c("func", "label", "shape")
      for (gg in (gii_ext_supported)) {
        if (endsWith(BOLD, paste0(".", gg, ".gii"))) { Bformat2 <- gg }
      }
      if (!is.na(Bformat2)) {
        Bformat2 <- switch(Bformat2, func="metric", label="label", shape="metric")
      } else {
        if (endsWith(BOLD, ".surf.gii")) { 
          Bformat2 <- "surf"
        }
      }
    } else if (any(endsWith(BOLD, c(".nii", ".nii.gz")))) {
      Bformat <- "NIFTI"
    }

  # xifti, gifti, or nifti data object
  } else if (inherits(BOLD, "xifti")) {
    Bformat <- "xifti"
    if (!is.null(BOLD$meta$cifti$intent)) {
      Bformat2 <- ciftiTools::supported_intents()$extension[match(
        BOLD$meta$cifti$intent,
        ciftiTools::supported_intents()$value
      )]
      Bformat2 <- gsub("\\.nii$", "", Bformat2)
    }
  } else if (inherits(BOLD, "gifti")) {
    Bformat <- "gifti"
    if(length(unique(BOLD$data_info$Intent))==1) {
      Bformat2 <- BOLD$data_info$Intent[1]
      if (length(unique(Bformat2))==1) {
        Bformat2 <- ifelse(Bformat2[1]=="NIFTI_INTENT_LABEL", "label", "metric")
      } else {
        if (length(Bformat2)==2 && all(Bformat2 == c("NIFTI_INTENT_POINTSET", "NIFTI_INTENT_TRIANGLE"))) {
          Bforamt2 <- "surf"
        }
      }
    }
  } else if (any(inherits(BOLD, nii_classes))) {
    Bformat <- "nifti"
    for (nn in nii_classes) {
      if (inherits(BOLD, nn)) { Bformat2 <- nn }
    }

  # List: GIFTI right and left
  } else if (is.list(BOLD) && (length(BOLD)==2)) {
    if (is.character(BOLD[[1]]) && is.character(BOLD[[2]])) {
      if (length(BOLD[[1]])==1 && length(BOLD[[2]])==1) {
        if (all(endsWith(c(BOLD[[1]], BOLD[[2]]), "gii"))) {
          Bformat <- "GIFTI2"
        }
      }
    } else if (inherits(BOLD[[1]], "gifti") && inherits(BOLD[[2]], "gifti")) {
      Bformat <- "gifti2"
    }

  # Data matrix
  } else if (is.numeric(BOLD)) {
    BOLD_dims_lens <- length(dim(BOLD))
    if (BOLD_dims_lens==4) {
      Bformat <- "nifti" # 4D array: treat as a "nifti"
    } else if (BOLD_dims_lens==2) {
      Bformat <- "data"
    } else if (BOLD_dims_lens==1) {
      warning("`BOLD` should not be vectorized.")
    } else {
      warning("`BOLD` should be TxV or a 4D array.")
    }
  } 
  
  # Print result.
  if (!is.na(Bformat)) {
    if (verbose) { cat("Inferred input format:", Bformat, "\n") }
  } else {
    warning("Could not infer BOLD format.")
  }

  # Return result.
  c(Bformat, Bformat2)
}