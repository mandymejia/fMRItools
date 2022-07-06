#' Convert CIFTI, NIFTI, or GIFTI file input to \eqn{T} by \eqn{V} matrix
#' 
#' Convert input data to \eqn{T} by \eqn{V} matrix. Reads in a CIFTI, NIFTI, 
#'  or GIFTI file. Transposes CIFTI and xifti object input. Masks NIFTI and 
#'  nifti object input.
#' 
#' @param x The object to coerce to a matrix
#' @param sortSub Sort subcortex by labels? Default: \code{FALSE}
#' @param verbose Print updates? Default: \code{FALSE}
#' @return x as a matrix.
#' @export
as.matrix_ifti <- function(x, sortSub=FALSE, verbose=FALSE) {

  format <- infer_format_ifti(x)
  if (format %in% c("CIFTI", "xifti")) {
    if (format == "CIFTI") {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read input data. Please install it", call. = FALSE)
      }
      x <- ciftiTools::read_cifti(x, brainstructures=ciftiTools::info_cifti(x)$cifti$brainstructures)
    }
    if (sortSub && !is.null(x$data$subcort)) {
      x$data$subcort <- x$data$subcort[order(x$meta$subcort$labels),]
    }
    x <- t(as.matrix(x))
  } else if (format %in% c("GIFTI", "gifti")) {
    if (format == "GIFTI") {
      if (!requireNamespace("gifti", quietly = TRUE)) {
        stop("Package \"gifti\" needed to read `X`. Please install it", call. = FALSE)
      }
      x <- gifti::read_gifti(x)
    }
    x <- t(do.call(cbind, x$data))
    if (verbose) {cat("GIFTI dimensions:\n"); print(dim(x))}
  } else if (format %in% c("NIFTI", "nifti")) {
    if (format == "NIFTI") {
      x <- read_nifti(x)
    }
    if (verbose) {cat("Masking NIFTI by removing locations with constant zero, NA, or NaN.\n")}
    z <- array(x %in% c(0, NA, NaN), dim=dim(x))
    mask <- !apply(z, seq(3), all)
    x <- matrix(x[rep(mask, dim(x)[4])], ncol=dim(x)[4])
    x <- t(x)
    if (verbose) {cat("NIFTI dimensions:\n"); print(dim(x))}
  }

  x
}