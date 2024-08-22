#' Generate cosine bases for the DCT
#'
#' https://en.wikipedia.org/wiki/Discrete_cosine_transform "DCT II"
#'
#' @param T_ Length of timeseries
#' @param n Number of cosine bases
#'
#' @return Matrix with cosine bases along columns
#'
#' @export
dct_bases <- function(T_, n){

  stopifnot(is_posNum(T_))
  stopifnot(is_posNum(n))
  stopifnot(n <= T_)

  b <- matrix(NA, T_, n)
  idx <- (seq(0, T_-1) + (1/2)) / (T_) # used to be `(seq(0, T_-1)) / (T_-1)`
  for (ii in seq(n)) {
    b[,ii] <- if (ii < T_) {
      cos(idx*pi*ii)
    } else {
      rep(c(2/pi, -2/pi), ceiling(T_/2))[seq(T_)] # cosine formula seems to break for this case.
    }
  }
  b
}

#' DCT and frequency conversion
#'
#' Convert between number of DCT bases and Hz of highpass filter
#'
#' Provide either \code{n} or \code{f} to calculate the other.
#'
#' If only the total length of the scan is known, you can set that to \code{TR}
#' and use \code{T_=1}.
#'
#' \eqn{f = n / (2 * T_ * TR)}
#'
#' @param T_ Length of timeseries (number of timepoints)
#' @param TR TR of the fMRI scan, in seconds (the time between timepoints)
#' @param n Number of cosine bases
#' @param f Hz of highpass filter
#'
#' @return If \code{n} was provided, the highpass filter cutoff (Hz) is returned.
#'  Otherwise, if \code{f} was provided, the number of cosine bases is returned.
#'  The result should be rounded before passing to \code{\link{dct_bases}}
#'
#' @export
dct_convert <- function(T_, TR, n=NULL, f=NULL){

  stopifnot(is_posNum(T_))
  stopifnot(is_posNum(TR))
  stopifnot(xor(is.null(n), is.null(f)))

  if (is.null(n)) {
    return(f * 2 * T_ * TR)
  } else if (is.null(f)) {
    return(n / (2 * T_ * TR))
  } else { stop() }
}

#' @rdname dct_convert
dct2Hz <- function(T_, TR, n){
  dct_convert(T_, TR, n=n)
}

#' @rdname dct_convert
Hz2dct <- function(T_, TR, f){
  dct_convert(T_, TR, f=f)
}

#' FFT detrending
#'
#' @param X \eqn{T \times V} numeric matrix. Each column is a voxel or vertex
#'  time series.
#' @param N Number of FFT entries to remove from each end of the vector
#'
#' @return Detrended \code{X}
#'
#' @importFrom stats mvfft
#'
#' @keywords internal
fft_detrend <- function(X, N) {

  stopifnot(is.numeric(X) && (length(dim(X))==2))
  stopifnot(is_posNum(N))

  T_ <- nrow(X)
  Y <- mvfft(X)
  Y[seq(N),] <- 0
  Y[seq(T_-N+1, T_),] <- 0
  Re(stats::mvfft(Y, inverse=TRUE)) / T_
}

#' Temporal filtering (bandpass, highpass, lowpass) with DCT or FFT
#'
#' @param X A numeric matrix, with each column being a timeseries to filter
#'  For fMRI data, \code{X} should be \code{T} timepoints by \code{V} brain
#'  locations.
#' 
#'  Alternatively, a single integer giving the number of timepoints in data.
#'  The return value will be the suitable set of DCT bases. Only works with
#'  \code{method == "DCT"}.
#' @param TR The time step between adjacent rows of \code{X}, in seconds.
#' @param hpf The frequency of the highpass filter, in Hertz. Default: \code{.008}.
#' @param lpf The frequency of the lowpass filter, in Hertz. Default: \code{NULL}
#'  (skip lowpass filtering). If both are provided, \code{lpf > hpf} must be true.
#' @param method \code{"DCT"} (default) or \code{"FFT"}. FFT is not compatible
#'  with \codE{lpf} yet.
#' @param verbose Print messages? Default: \code{FALSE}.
#'
#' @return Filtered \code{X}, or if \code{X} was an integer, the set of DCT
#'  bases to use for nuisance regression (not including an intercept).
#'
#' @export
#'
#' @examples
#' temporal_filter(matrix(rnorm(700), nrow=100), TR=.72)
temporal_filter <- function(
  X, TR,
  hpf=.008, lpf=NULL,
  method=c("DCT", "FFT"),
  verbose=FALSE) {

  # Argument checks.
  X_int <- is_posNum(X) && X == round(X)
  if (!X_int) {
    X <- as.matrix(X)
    stopifnot(is.numeric(X))
  }
  stopifnot(is_posNum(TR))
  stopifnot(is.null(hpf) || is_posNum(hpf))
  stopifnot(is.null(lpf) || is_posNum(lpf))
  if (!is.null(hpf)) { stopifnot(hpf < 1/TR/2) }
  if (!is.null(hpf) && !is.null(lpf)) { stopifnot(lpf > hpf) }
  if (is.null(hpf) && is.null(lpf)) { return(colCenter(X)) }
  method <- match.arg(method, c("DCT", "FFT"))

  # Get number of timepoints.
  T_ <- if (X_int) { X } else { nrow(X) }

  # DCT ------------------------------------------------------------------------
  if (method == "DCT") {
    # Get vector of frequencies associated with each DCT basis.
    # 0 and T+1 values used to add opt to skip filtering, or raising error.
    freq_exp <- dct_convert(T_, TR, n=seq(0, T_+1)) 
    # Compute the bases. 
    bases <- dct_bases(T_, T_)

    # Initialize nuisance matrix.
    nmat <- NULL

    # Collect DCT HPF bases.
    if (!is.null(hpf)) {
      hpf_idx <- which.min(abs(freq_exp - hpf)) - 1 # 0 or T_: skip/error
      if (hpf_idx == 0) {
        if (verbose) { cat("0 bases for HPF (after rounding).\n") }
      } else if (hpf_idx == T_+1) {
        stop("`hpf` is too high.")
      } else {
        if (verbose) { cat(hpf_idx, "bases for HPF.\n") }
        if (hpf_idx == T_) { warning("No DOF left.") }
        nmat_hpf <- bases[,seq(hpf_idx),drop=FALSE]
        colnames(nmat_hpf) <- paste(round(freq_exp[1+seq(hpf_idx)], 4), "Hz")
        nmat <- cbind(nmat, nmat_hpf)
      }
    }

    # Collect the DCT LPF bases.
    if (!is.null(lpf)) {
      lpf_idx <- which.min(abs(freq_exp - lpf)) - 1 # 0 or T_: error/skip
      # Check redundancy with hpf.
      if (!is.null(hpf)) {
        if (lpf_idx == hpf_idx) {
          lpf_idx <- lpf_idx + 1 # avoid regressing same idx twice.
        }
        if (lpf_idx == hpf_idx + 1) {
          warning("No DOF after LPF and HPF.")
        }
      }
      if (lpf_idx == 0) {
        stop("`lpf` is too low.")
      } else if (lpf_idx == T_+1) {
        if (verbose) { cat("0 bases for LPF (after rounding).\n") }
      } else {
        if (verbose) { cat(T_ - lpf_idx + 1, "bases for LPF.\n") }
        if (lpf_idx == 1) { warning("No DOF left.") }
        nmat_lpf <- bases[,seq(lpf_idx, T_),drop=FALSE]
        colnames(nmat_lpf) <- paste(round(freq_exp[1+seq(lpf_idx, T_)], 4), "Hz")
        nmat <- cbind(nmat, nmat_lpf)
      }
    }

    if (X_int) { return(nmat) } 

    # Do nuisance regression with an intercept and the collected DCT bases.
    nmat <- cbind(1, nmat)
    out <- nuisance_regression(X, nmat)

  # FFT ------------------------------------------------------------------------
  } else if (method == "FFT") {
    if (X_int) { stop("Not applicable: integer X value for FFT (number of timepoints). Please provide data.") }
    if (!is.null(lpf)) { stop("Not implemented: lowpass-filtering with FFT.") }

    # [NOTE] This code used to be used for DCT too, before lpf was implemented.
    N <- switch(method,
      #DCT = round(hpf * TR * T_ * 2),
      FFT = round(hpf * TR * T_)
    )
    if (N < 1) { return(colCenter(X)) }
    if (N > T_/2) { #ifelse(method=="DCT", T_, T_/2)) {
      stop("Maximum `f` for this data is: ", round(1/(TR*2), digits=5), " Hz.")
    }
    out <- switch(method,
      #DCT = nuisance_regression(X, cbind(1, dct_bases(T_, N))),
      FFT = fft_detrend(X, N)
    )
  }

  out
}

#' @rdname temporal_filter
#' @export
detrend <- function(X, TR, hpf=.008, method=c("DCT", "FFT")) {
  temporal_filter(X=X, TR=TR, hpf=hpf, lpf=NULL, method=method)
}