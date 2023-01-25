#' Carpetplot
#' 
#' Plot a matrix with \code{graphics::image}. For fMRI data, this is the 
#'  "carpetplot" or grayplot coined by (Power, 2017). The \code{graphics} and
#'  \code{grDevices} packages are required.
#' 
#' @param x The \eqn{T \times V} numeric data matrix, or a \code{"xifti"} object.
#'  In the plot, the \eqn{T} index will increase from left to right, and the 
#'  \eqn{V} will increase from top to bottom.
#' @param qcut Sets blackpoint at the \code{qcut} quantile, and the
#'  whitepoint at the \code{1-qcut} quantile. Default: \code{.1}. This is
#'  equivalent to setting the color range between the 10% and 90% quantiles.
#'  The quantiles are computed across the entire data matrix after any 
#'  centering or scaling. 
#' 
#'  Must be between 0 and .49. If \code{0} or \code{NULL} (default), do not 
#'  clamp the data values.
#' @param fname A \code{.pdf} (highly recommended) or \code{.png} file path
#'  to write the carpetplot to. If \code{NULL} (default), return the plot directly
#'  instead of writing a file.
#' @param center,scale Center and scale the data? If \code{x} is fMRI data 
#'  which has not otherwise been centered or scaled, it is recommended to center
#'  but not scale it (default).
#' @param colors \code{"gray255"} (default) will use a grayscale color ramp
#'  from black to white. Otherwise, this should be a character vector of
#'  color names to use. 
#' 
#'  Colors will be assigned from the lowest to the highest data value, after 
#'  any clamping of the data values by \code{qcut}.
#' @param sortSub If \code{x} is a \code{"xifti"} object with subcortical data,
#'  should the voxels be sorted by structure alphabetically? Default: \code{TRUE}.
#' @param ... Additional arguments to \code{pdf} or \code{png}, such as width
#'  and height.
#' 
# @importFrom grDevices dev.off pdf png gray.colors
#' 
#' @importFrom stats quantile
#' @return The image or \code{NULL, invisibly} if a file was written.
#' 
#' @section References:
#'  \itemize{
#'    \item{Power, J. D. A simple but useful way to assess fMRI scan qualities. NeuroImage 154, 150-158 (2017).}
#' }
#' 
#' @export
#' 
carpetplot <- function(
  x, qcut=.1, fname=NULL, center=TRUE, scale=FALSE, colors="gray255", sortSub=TRUE, ...){

  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed since `svd` failed. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package \"graphics\" needed since `svd` failed. Please install it.", call. = FALSE)
  }

  # Get T x V matrix.
  x <- as.matrix_ifti(x, sortSub=sortSub)

  # Center or scale.
  if (center | scale) { x <- scale(x, center=center, scale=scale) }

  # Orient correctly for `image`
  x <- x[,seq(ncol(x), 1)]

  # Clamp to quantile cutoffs.
  if (!is.null(qcut)) {
    qcut <- as.numeric(qcut)
    if (qcut > 1e-8) {
      stopifnot(qcut < .49)
      x[x < quantile(x, qcut)] <- quantile(x, qcut)
      x[x > quantile(x, 1-qcut)] <- quantile(x, 1-qcut)
    }
  }

  # Check colors.
  if (identical(colors, "gray255")) { 
    colors <- grDevices::gray.colors(255, start=0, end=1)
  }

  # Check file name
  if (!is.null(fname)) {
    fname <- as.character(fname)
    ftype <- ifelse(endsWith(fname, ".png"), "png", "pdf")
    if (ftype == "pdf" && !endsWith(fname, ".pdf")) {
      fname <- paste0(fname, ".pdf")
    }

    switch(ftype, png=grDevices::png, pdf=grDevices::pdf)(fname, ...)
  }

  graphics::image(
    x, 
    axes=FALSE, frame.plot = FALSE, xlab="", ylab="", ann=FALSE, 
    useRaster=TRUE, col=colors
  )

  if (!is.null(fname)) { 
    grDevices::dev.off()
  }

  invisible(NULL)
}

#' Stacked carpetplot
#' 
#' Stacks carpetplots on top of one another by rbinding the matrices. 
#' 
#' @param x_list List of data matrices
#' @param center,scale Center and scale the data? If \code{x} is fMRI data 
#'  which has not otherwise been centered or scaled, it is recommended to center
#'  but not scale it (default).
#' @param qcut Sets blackpoint at the \code{qcut} quantile, and the
#'  whitepoint at the \code{1-qcut} quantile. Default: \code{.1}. This is
#'  equivalent to setting the color range between the 10% and 90% quantiles.
#'  The quantiles are computed across the entire data matrix after any 
#'  centering or scaling. 
#' 
#'  Must be between 0 and .49. If \code{0} or \code{NULL} (default), do not 
#'  clamp the data values.
#' @param match_scale Match the scales of the carpetplots? Default: \code{TRUE}.
#' @param nsep Equivalent number of data locations for size of gap between
#'  carpetplots. Default: zero (no gap).
#' @param ... Additional arguments to \code{carpetplot}
#' @return \code{NULL}, invisibly
#' @export
carpetplot_stack <- function(x_list, center=TRUE, scale=FALSE, qcut=.1, match_scale=TRUE, nsep=0, ...){

  # [TO DO] check args
  args <- list(...)
  sortSub <- ifelse("sortSub" %in% names(args), args$sortSub, TRUE)

  x_list <- lapply(x_list, as.matrix_ifti, sortSub=sortSub)

  # Center or scale.
  if (center | scale) { 
    x_list <- lapply(x_list, scale, center=center, scale=scale)
  }

  # Clamp to quantile cutoffs.
  if (!is.null(qcut)) {
    qcut <- as.numeric(qcut)
    if (qcut > 1e-8) {
      stopifnot(qcut < .49)
      for (ii in seq(length(x_list))) {
        x_list[[ii]][x_list[[ii]] < quantile(x_list[[ii]], qcut)] <- quantile(x_list[[ii]], qcut)
        x_list[[ii]][x_list[[ii]] > quantile(x_list[[ii]], 1-qcut)] <- quantile(x_list[[ii]], 1-qcut)
      }
    }
  }

  if (match_scale) {
    scale_cp <- function(x){ q <- min(x[]); x[] <- (x[]-q)/(max(x[]) - q); x }
    x_list <- lapply(x_list, scale_cp)
  }
  xmax <- max(vapply(x_list, max, 0, na.rm=TRUE), na.rm=TRUE)

  if (nsep > 0) {
    x2 <- vector("list", length(x_list)*2 - 1)
    for (ii in seq(length(x2))) {
      if (ii %% 2 == 1) {
        x2[[ii]] <- x_list[[ceiling(ii/2)]]
      } else {
        x2[[ii]] <- matrix(xmax, nrow=nrow(x_list[[1]]), ncol=nsep)
      }
    }
    x2 <- do.call(cbind, x2)
  } else {
    x2 <- do.call(cbind, x_list)
  }

  carpetplot(x2, center=FALSE, scale=FALSE, qcut=0, ...)
}