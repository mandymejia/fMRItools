#' \code{bptf} function from FSL
#'
#' Copy of \code{bptf} highpass filter from FSL. The results are very similar
#'  but not identical.
#'
#' Sources:
#'  https://cpb-us-w2.wpmucdn.com/sites.udel.edu/dist/7/4542/files/2016/09/fsl_temporal_filt-15sywxn.m
#'  https://github.com/rordenlab/niimath/blob/master/src/coreFLT.c#L1935
#'
#' @param orig_data \eqn{T \times V} data matrix whose columns will be detrended
#' @param HP_sigma The frequency parameter, sigma, for the highpass filter
#' @param LP_sigma The frequency parameter, sigma, for the lowpass filter
#'
#' @return The data with detrended columns
#'
#' @export
#'
#' @section References:
#'  \itemize{
#'    \item{Jenkinson, M., Beckmann, C. F., Behrens, T. E. J., Woolrich, M. W. & Smith, S. M. FSL. NeuroImage 62, 782-790 (2012).}
#' }
#'
#' @examples
#' fsl_bptf(matrix(rnorm(700), nrow=100))
fsl_bptf <- function(orig_data, HP_sigma=2000, LP_sigma=NULL) {

  orig_data <- as.matrix(orig_data)
  stopifnot(is.numeric(orig_data))
  nT <- nrow(orig_data)

  # Coerce `HP_sigma` and `LP_sigma` to NULL or a positive number.
  if (!is.null(HP_sigma)) {
    stopifnot(is_1(HP_sigma))
    if (HP_sigma<=0) { HP_sigma <- NULL }
  }
  if (!is.null(LP_sigma)) {
    stopifnot(is_1(LP_sigma))
    if (LP_sigma<=0) { LP_sigma <- NULL }
  }
  do_HPF <- !is.null(HP_sigma)
  do_LPF <- !is.null(LP_sigma)

  if (do_HPF) {
    orig_data <- nuisance_regression(orig_data, cbind(1, seq(nT)))
    
    HP_filt_size <- ceiling(HP_sigma*3)#round(HP_sigma*8)
    HP_lin <- seq(-HP_filt_size/2, HP_filt_size/2, length.out=HP_filt_size)
    HP_gfilt <- exp( -(HP_lin^2) / (2*(HP_sigma^2)) )
    HP_gfilt <- HP_gfilt/sum(HP_gfilt)
  }

  if (do_LPF) {
    LP_filt_size <- ceiling(LP_sigma*3)#round(LP_sigma*8)
    LP_lin <- seq(-LP_filt_size/2, LP_filt_size/2, length.out=LP_filt_size)
    LP_gfilt <- exp( -(LP_lin^2) / (2*(LP_sigma^2)) )
    LP_gfilt <- LP_gfilt/sum(LP_gfilt)
  }

  if (do_HPF) {
    filt_data <- matrix(NA, nrow=nT, ncol=ncol(orig_data))
    back <- floor((HP_filt_size-1)/2)
    front <- ceiling((HP_filt_size-1)/2)
    for (t in seq(nT)) {
      if ((t-back < 1) && (t+front > nT)) {
        trunc_HP_gfilt <- HP_gfilt[seq(back-t+2, HP_filt_size-(t+front-nT))]
        trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
        filt_data[t,] <- trunc_HP_gfilt %*% orig_data
      } else if (t-back < 1) {
        trunc_HP_gfilt <- HP_gfilt[seq(back-t+2, HP_filt_size)]
        trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
        filt_data[t,] <- trunc_HP_gfilt %*% orig_data[seq(t+front),,drop=FALSE]
      } else if (t+front > nT) {
        trunc_HP_gfilt <- HP_gfilt[seq(HP_filt_size-(t+front-nT))]
        trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
        filt_data[t,] <- trunc_HP_gfilt %*% orig_data[seq(t-back, nT),,drop=FALSE]
      } else {
        filt_data[t,] <- HP_gfilt %*% orig_data[seq(t-back, t+front),,drop=FALSE]
      }
    }
    filt_data <- orig_data - filt_data
  } else {
    filt_data <- orig_data
  }

  if (do_LPF) {
    filt_data_orig <- filt_data
    back <- floor((LP_filt_size-1)/2)
    front <- ceiling((LP_filt_size-1)/2)
    for (t in seq(nT)) {
      if ((t-back < 1) && (t+front > nT)) {
        trunc_LP_gfilt <- LP_gfilt[seq(back-t+2, LP_filt_size-(t+front-nT))]
        trunc_LP_gfilt <- trunc_LP_gfilt/sum(trunc_LP_gfilt)
        filt_data[t,] <- trunc_LP_gfilt %*% filt_data_orig[seq(t-back, t+front),,drop=FALSE]
      } else if (t-back < 1) {
        trunc_LP_gfilt <- LP_gfilt[seq(back-t+2, LP_filt_size)]
        trunc_LP_gfilt <- trunc_LP_gfilt/sum(trunc_LP_gfilt)
        filt_data[t,] <- trunc_LP_gfilt %*% filt_data_orig[seq(t+front),,drop=FALSE]
      } else if (t+front > nT) {
        trunc_LP_gfilt <- LP_gfilt[seq(LP_filt_size-(t+front-nT))]
        trunc_LP_gfilt <- trunc_LP_gfilt/sum(trunc_LP_gfilt)
        filt_data[t,] <- trunc_LP_gfilt %*% filt_data_orig[seq(t-back, nT),,drop=FALSE]
      } else {
        filt_data[t,] <- LP_gfilt %*% filt_data_orig[seq(t-back, t+front),,drop=FALSE]
      }
    }
  }

  filt_data
}