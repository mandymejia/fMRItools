#' Convert CIFTI, NIFTI, or GIFTI input to \eqn{T \times V} matrix
#' 
#' Convert CIFTI, NIFTI, or GIFTI input to a \eqn{T \times V} matrix by 
#'  reading it in with the corresponding package, and then separating the data
#'  from the metadata. Also works with the intermediate R objects created from
#'  reading these files: \code{"xifti"} objects from \code{ciftiTools}, 
#'  \code{"gifti"} objects from \code{gifti}, 
#'  \code{"nifti"} or \code{"niftiExtension"} objects from \code{oro.nifti}, and
#'  \code{"niftiImage"} objects from \code{RNifti}.
#'  
#'  For CIFTI files, only intents supported by \code{ciftiTools} are supported: 
#'  \code{dscalar}, \code{dtseries}, and \code{dlabel}. For NIFTI file or
#'  NIFTI-intermediate R object input, the data will be vectorized/masked. 
#' 
#' @param x The object to coerce to a matrix
#' @param meta Return metadata too? Default: \code{FALSE}.
#' @param sortSub For CIFTI format input only. Sort subcortex by labels? 
#'  Default: \code{FALSE}
#' @param TbyV Return the data matrix in \eqn{T \times V} form? Default:
#'  \code{TRUE}. If \code{FALSE}, return in \eqn{V \times T} form instead.
#'  Using this argument may be faster than transposing after the function call. 
#' @param verbose Print updates? Default: \code{FALSE}
#' @param ... If \code{x} is a file path, additional arguments to the function
#'  used to read in \code{x} can be specified here. For example, if \code{x}
#'  is a path to a CIFTI file, \code{...} might specify which \code{idx} and
#'  \code{brainstructures} to read in.
#' 
#' @return If \code{!meta}, \code{x} as a matrix. If \code{meta}, a list of
#'  length two: the first entry is \code{x} as a matrix, and the second entry is 
#'  the metadata of \code{x}. 
#' 
#' @export
as.matrix_ifti <- function(
  x, meta=FALSE, sortSub=FALSE, TbyV=TRUE, verbose=FALSE, ...) {

  x_meta <- NULL

  stopifnot(is_1(meta, "logical"))
  stopifnot(is_1(sortSub, "logical"))
  stopifnot(is_1(TbyV, "logical"))
  stopifnot(is_1(verbose, "logical"))

  # Get the format of `x`
  format <- infer_format_ifti(x)
  format2 <- format[2]; format <- format[1]

  # Handle CIFTI input.
  if (format %in% c("CIFTI", "xifti")) {
    if (format == "CIFTI") {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read input data. Please install it", call. = FALSE)
      }
      x_readArgs <- list(...)
      if (!("brainstructures" %in% names(x_readArgs))) {
        x_readArgs$brainstructures <- ciftiTools::info_cifti(x)$cifti$brainstructures
      }
      x <- do.call(ciftiTools::read_cifti, c(list(x), x_readArgs))
    }
    if (sortSub && !is.null(x$data$subcort)) {
      x$data$subcort <- x$data$subcort[order(x$meta$subcort$labels),]
    }
    if (meta) {
      x_meta <- c(
        list(brainstructures_nV=lapply(x$data, nrow)), 
        x[names(x)[names(x)!="data"]]
      )
    }
    x <- as.matrix(x)
    if (TbyV) { x <- t(x) }

  # Handle GIFTI input.
  } else if (format %in% c("GIFTI", "gifti")) {
    if (!is.na(format2) && format2 == "surf") { stop(
      "`x` represents surface data, not BOLD/fMRI data."
    ) }
    if (format == "GIFTI") {
      if (!requireNamespace("gifti", quietly = TRUE)) {
        stop("Package \"gifti\" needed to read `X`. Please install it", call. = FALSE)
      }
      x <- gifti::read_gifti(x)
    }
    if (meta) {
      x_meta <- x[names(x)[names(x)!="data"]]
    }
    x <- t(do.call(cbind, x$data))
    if (verbose) {cat("GIFTI dimensions:\n"); print(dim(x))}

  # Handle NIFTI input.
  } else if (format %in% c("NIFTI", "nifti")) {
    if (format == "NIFTI") {
      x <- read_nifti(x)
    }
    if (verbose) {cat("Masking NIFTI by removing locations with constant zero, NA, or NaN.\n")}
    z <- array(x %in% c(0, NA, NaN), dim=dim(x))
    mask <- !apply(z, seq(3), all)
    if (meta) {
      x_meta <- list(mask=mask, attributes=attributes(x))
    }
    x <- matrix(x[rep(mask, dim(x)[4])], ncol=dim(x)[4])
    x <- t(x)
    if (verbose) {cat("NIFTI dimensions:\n"); print(dim(x))}
  }

  # CIFTI: add sortSub metadata
  out <- if (meta) {
    list(data=x, meta=x_meta)
  } else {
    x
  }
}