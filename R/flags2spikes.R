#' Flags to nuisance spikes
#' 
#' Convert flagged volumes to corresponding one-hot encoded vectors
#'  which can be used for nuisance regression.
#' 
#' @param flags Numeric vector of integers indicating the indices of 
#'  the flagged volumes. Or, a logical vector of length \code{n_time} where 
#'  \code{TRUE} values indicate the flagged volumes.
#' @param n_time The length of the vectors to obtain. For nuisance regression,
#'  this is the length of the BOLD data. The highest index in \code{flags}
#'  should not exceed \code{n_time}.
#' 
#' @return A numeric matrix of ones and zeroes. The number of rows will be
#'  \code{n_time} and the number of columns will be the number of flags. Each 
#'  column will have a \code{1} at the flag index, and \code{0} elsewhere.
#' @export
flags2spikes <- function(flags, n_time){
  stopifnot(is_posNum(n_time) && n_time==round(n_time))
  if (is.logical(flags)) {
    stopifnot(length(flags) == n_time, !anyNA(flags))
    flags <- which(flags)
  }
  stopifnot(is.numeric(flags), !anyNA(flags))
  stopifnot(all(flags == round(flags)))

  flags <- sort(unique(flags))
  nF <- length(flags)
  if (nF == 0) { return(matrix(0, nrow=n_time, ncol=0)) }
  stopifnot(flags[1] >= 1 && flags[nF] <= n_time)

  spikes <- matrix(0, nrow=n_time, ncol=nF)
  spikes[cbind(flags, seq_len(nF))] <- 1
  spikes
}
