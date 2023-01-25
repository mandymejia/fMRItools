#' Bandstop filter
#' 
#' Filter out frequencies within a given range using a Chebyshev Type II 
#'  stopband. Essentially a convenience wrapper for the 
#'  \code{\link[gsignal]{cheby2}} function. 
#' 
#' @param X A numeric matrix, with each column being a timeseries to apply
#'  the stopband filter. For fMRI data, \code{X} should be \code{T} timepoints
#'  by \code{V} brain locations.
#' @param TR The time step between adjacent rows of \code{x}, in seconds
#' @param f1,f2 The frequency limits for the filter, in Hz. \code{f1 < f2}
#' @param Rs The amount of attenuation of the stopband ripple, in dB
#' 
#' @return The filtered data
#' @export
#' 
#' @examples
#' if (requireNamespace("gsignal", quietly=TRUE)) {
#'  n_voxels = 1e4
#'  n_timepoints = 100
#'  X = cbind(arima.sim(n=100, list(ar=.6)), arima.sim(n=100, list(ar=.6)))
#'  Y = bandstop_filter(X, .72, .31, .43)
#' }
bandstop_filter <- function(X, TR, f1, f2, Rs=20){

  if (!requireNamespace("gsignal", quietly = TRUE)) {
    stop("Package \"gsignal\" needed for bandstop filter.", call. = FALSE)
  }

  if (is.vector(X)) { X <- as.matrix(X) }
  stopifnot(nrow(X)>5)
  stopifnot(is_posNum(TR))
  stopifnot(is_posNum(f1))
  stopifnot(is_posNum(f2))
  stopifnot(f1 < f2)

  nf <- c(f1, f2) * TR * 2
  nfilt <- gsignal::cheby2(2, Rs=Rs, w=nf, type="stop")
  gsignal::filtfilt(nfilt, X)
}