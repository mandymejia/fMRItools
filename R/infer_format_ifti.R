#' Infer fMRI data format
#'
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return A length-two vector. The first element indicates the format:
#'  \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"GIFTI"} file path, \code{"gifti"} object,
#'  \code{"NIFTI"} file path, \code{"nifti"} object,
#'  \code{"RDS"} file path, or \code{"data"}.
#'  The second element indicates
#'  the sub-format if relevant; i.e. the type of CIFTI or GIFTI file/object.
#'
#' @export
infer_format_ifti <- function(BOLD, verbose=FALSE){

  Bformat <- Bformat2 <- NA

  # Define some file extensions to look out for
  # CIFTI: copied from `ciftiTools::supported_intents()`
  cii_supp_int <- data.frame(rbind(
    c("dtseries.nii", "NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES",  3002, "ConnDenseSeries"),
    c("dscalar.nii",  "NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS", 3006, "ConnDenseScalar"),
    c("dlabel.nii",   "NIFTI_INTENT_CONNECTIVITY_DENSE_LABELS",  3007, "ConnDenseLabel")
  ))
  colnames(cii_supp_int) <- c("extension", "intent_code", "value", "intent_name")
  cii_ext_supported <- paste0(".", cii_supp_int$extension)
  cii_ext_notsupport <- paste0(
    ".", c(
      "dconn", "pconn", "ptseries", "dtraj", "pscalar",
      "pdconn", "dpconn", "pconnseries", "pconnscalar"
    ), ".nii"
  )
  # GIFTI
  gii_ext_supported <- c("func", "label")
  # NIFTI
  nii_classes <- c(oro.nifti="niftiExtension", RNifti="niftiImage", oro.nifti="nifti")

  # Character: CIFTI, GIFTI, NIFTI, or RDS
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
    } else if (any(endsWith(BOLD, c(".rds", ".RDS")))) {
      Bformat <- "RDS"
    }

  # xifti, gifti, or nifti data object
  } else if (inherits(BOLD, "xifti")) {
    Bformat <- "xifti"
    if (!is.null(BOLD$meta$cifti$intent)) {
      Bformat2 <- cii_supp_int$extension[match(
        BOLD$meta$cifti$intent,
        cii_supp_int$value
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
      warning("`BOLD` should not be vectorized, if it's numeric data.")
    } else {
      warning("`BOLD` should be TxV or a 4D array, if it's numeric data.")
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

#' Infer fMRI data format for several inputs
#'
#' Vectorized version of \code{\link{infer_format_ifti}}. Expects all inputs
#'  to have the same format.
#'
#' Raises an error if the elements of \code{BOLD} do not share the same format.
#'
#' @param BOLD The vector of fMRI data, expected to be of one format
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return A length-two vector. The first element indicates the format:
#'  \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"GIFTI"} file path, \code{"gifti"} object,
#'  \code{"NIFTI"} file path, \code{"nifti"} object,
#'  \code{"RDS"} file path, or \code{"data"}. The second element indicates
#'  the sub-format if relevant; i.e. the type of CIFTI or GIFTI file/object.
#'
#' @export
infer_format_ifti_vec <- function(BOLD, verbose=FALSE){
  BOLD <- as.list(BOLD)
  Bformat <- lapply(BOLD, infer_format_ifti)
  Bformat <- unique(Bformat)
  if (length(Bformat)>1) {
    stop(paste(
      "The formats are not identical: ",
      paste(lapply(Bformat, paste, collapse=", "), collapse="; ")
    ))
  }
  if (verbose) { cat("Inferred input format:", Bformat[[1]], "\n") }
  Bformat[[1]]
}
