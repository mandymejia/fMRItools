#' Bandstop filter
#' 
#' Filter out frequencies within a given range using a Chebyshev Type II 
#'  stopband. Basically a convenience wrapper for the \code{gsignal::cheby2}
#'  function. 
#' 
#' @param x A numeric matrix, with each column being a timeseries to apply
#'  the stopband filter
#' @param TR the time difference between rows in x, in seconds.
#' @param f1,f2 The frequency limits for the filter, in Hz. \code{f1 < f2}.
#' @param Rs The amount of attenuation of the stopband ripple, in dB
#' 
#' @return The filtered data
#' @export
#' 
#' @examples
#' library(gsignal)
#' n_voxels = 1e4
#' n_timepoints = 100
#' X = cbind(arima.sim(n=100, list(ar=.6)), arima.sim(n=100, list(ar=.6)))
#' Y = bandstop_filter(X, .72, .31, .43)
bandstop_filter <- function(x, TR, f1, f2, Rs=20){

  if (!requireNamespace("gsignal", quietly = TRUE)) {
    stop("Package \"gsignal\" needed for bandstop filter.", call. = FALSE)
  }

  if (is.vector(x)) { x <- as.matrix(x) }
  stopifnot(nrow(x)>5)
  is_pos_num <- function(q){ length(q)==1 && is.numeric(q) && q>0 }
  stopifnot(is_pos_num(TR))
  stopifnot(is_pos_num(f1)); stopifnot(is_pos_num(f2))
  stopifnot(f1 < f2)

  nf <- c(f1, f2) * TR * 2
  nfilt <- gsignal::cheby2(2, Rs=Rs, w=nf, type="stop")
  gsignal::filtfilt(nfilt, x)
}