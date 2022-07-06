#' \code{bptf} function from FSL
#' 
#' Copy of highpass filter as implemented by \code{bptf} in FSL. The results are very
#'  similar but not identical.
#' 
#' Sources:
#'  https://cpb-us-w2.wpmucdn.com/sites.udel.edu/dist/7/4542/files/2016/09/fsl_temporal_filt-15sywxn.m
#'  https://github.com/rordenlab/niimath/blob/master/src/core32.c
#' 
#' @param orig_data \eqn{T} by \eqn{V} data matrix whose columns will be detrended
#' @param HP_sigma The frequency parameter for the highpass filter
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
fsl_bptf <- function(orig_data, HP_sigma=2000) {
  T_ <- nrow(orig_data)

  orig_data <- nuisance_regression(orig_data, cbind(1, seq(T_)))

  HP_filt_size <- ceiling(HP_sigma*3)#round(HP_sigma*8)
  HP_lin <- seq(-HP_filt_size/2, HP_filt_size/2, length.out=HP_filt_size)
  HP_gfilt <- exp( -(HP_lin^2) / (2*(HP_sigma^2)) )
  HP_gfilt <- HP_gfilt/sum(HP_gfilt)

  filt_data <- matrix(NA, nrow=T_, ncol=ncol(orig_data))
  back <- floor((HP_filt_size-1)/2)
  front <- ceiling((HP_filt_size-1)/2)
  for (t in seq(T_)) {
    if ((t-back < 1) && (t+front > T_)) {
      trunc_HP_gfilt <- HP_gfilt[seq(back-t+2, HP_filt_size-(t+front-T_))]
      trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
      filt_data[t,] <- trunc_HP_gfilt %*% orig_data
    } else if (t-back < 1) {
      trunc_HP_gfilt <- HP_gfilt[seq(back-t+2, HP_filt_size)]
      trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
      filt_data[t,] <- trunc_HP_gfilt %*% orig_data[seq(t+front),]
    } else if (t+front > T_) {
      trunc_HP_gfilt <- HP_gfilt[seq(HP_filt_size-(t+front-T_))]
      trunc_HP_gfilt <- trunc_HP_gfilt/sum(trunc_HP_gfilt)
      filt_data[t,] <- trunc_HP_gfilt %*% orig_data[seq(t-back, T_),]
    } else {
      filt_data[t,] <- HP_gfilt %*% orig_data[seq(t-back, t+front),]
    }
  }
  orig_data - filt_data
}