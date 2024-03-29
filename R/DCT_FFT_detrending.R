#' Generate cosine bases for the DCT
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

  b <- matrix(NA, T_, n)
  idx <- (seq(T_)-1)/(T_-1)
  for (ii in seq(n)) { b[,ii] <- cos(idx*pi*ii) }
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

#' Detrending with DCT or FFT
#' 
#' @param X A numeric matrix, with each column being a timeseries to detrend. 
#'  For fMRI data, \code{X} should be \code{T} timepoints by \code{V} brain 
#'  locations.
#' @param TR The time step between adjacent rows of \code{X}, in seconds
#' @param f The frequency of the highpass filter, in Hertz. Default: \code{.008}
#' @param method \code{"DCT"} (default) or \code{"FFT"}.
#' 
#' @return Detrended \code{X}
#' 
#' @export 
#' 
#' @examples 
#' detrend(matrix(rnorm(700), nrow=100), TR=.72)
detrend <- function(X, TR, f=.008, method=c("DCT", "FFT")) {

  X <- as.matrix(X)
  stopifnot(is.numeric(X))
  stopifnot(is_posNum(TR))
  stopifnot(is_posNum(f))
  method <- match.arg(method, c("DCT", "FFT"))
  
  T_ <- nrow(X)
  N <- round(f * TR * T_)

  N <- switch(method,
    DCT = round(f * TR * T_ * 2),
    FFT = round(f * TR * T_)
  )

  if (N < 1) { return(t( scale(t(X), scale=FALSE)) ) }
  if (N > ifelse(method=="DCT", T_, T_/2)) { 
    stop("Maximum `f` for this data is: ", round(1/(TR*2), digits=5), " Hz.") 
  }

  switch(method,
    DCT = nuisance_regression(X, cbind(1, dct_bases(T_, N))),
    FFT = fft_detrend(X, N)
  )
}