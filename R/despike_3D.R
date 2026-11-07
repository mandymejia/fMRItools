#' 3dDespike from AFNI, step 0
#' 
#' Fourier basis for the smooth curve that \code{\link{despike_3D}} is based on.
#'  The number of basis functions is \code{2*round(nT/30/2)+1}.
#' 
#' @param nT The number of timepoints.
#' 
#' @return The \eqn{nT \times (2K+1)} basis matrix.
#' 
#' @keywords internal
despike_3D.basis <- function(nT){
  K <- round(nT/30/2)
  if (K < 1) { return(matrix(1, nT, 1)) }
  ang <- 2*pi*outer(seq_len(nT), seq_len(K))/nT
  cbind(1, sin(ang), cos(ang))
}

#' 3dDespike from AFNI, step 1
#' 
#' Compute the quantile regression that \code{\link{despike_3D}} is based on.
#' 
#' @param Yt The data vector.
#' @param basis Optional precomputed \code{\link{despike_3D.basis}}. Supply
#'  this when despiking many data vectors of equal length.
#' 
#' @return the quantile regression
#' 
#' @keywords internal
despike_3D.qreg <- function(Yt, basis=NULL){

  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package \"quantreg\" needed. Please install it", call. = FALSE)
  }

  if (is.null(basis)) { basis <- despike_3D.basis(length(Yt)) }
  stopifnot(nrow(basis) == length(Yt))
  quantile.reg <- quantreg::rq(Yt~basis-1)
}

#' 3dDespike from AFNI, step 2
#' 
#' Identify and interpolate the outliers for \code{\link{despike_3D}}.
#' 
#' @param qreg the quantile regression from \code{\link{despike_3D.qreg}}
#' @param c1 spike threshold. Default: \code{2.5}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. 
#' 
#' @importFrom stats residuals mad fitted
#' @keywords internal
despike_3D.interpolate <- function(qreg, c1=2.5, c2=4){
  qreg_resid <- residuals(qreg)
  mad_resid <- mad(qreg_resid)
  if (!is.finite(mad_resid) || mad_resid == 0) {
    return(fitted(qreg) + qreg_resid)
  }
  s <- qreg_resid/mad_resid
  a <- abs(s)
  spike <- a > c1
  a[spike] <- c1 + (c2-c1)*tanh((a[spike]-c1)/(c2-c1))
  fitted(qreg) + sign(s)*a*mad_resid
}

#' 3dDespike from AFNI
#' 
#' Identify and interpolate outliers. See 
#'  \href{https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html}{the
#'  AFNI documentation for 3dDespike} for additional information.
#' 
#' @param Yt The data vector.
#' @param c1 spike threshold. Default: \code{2.5}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. 
#' @param basis Optional precomputed \code{\link{despike_3D.basis}}.
#' 
#' @examples
#' if (requireNamespace("quantreg", quietly=TRUE)) {
#'  y <- rnorm(99) + cos(seq(99)/15)*3
#'  y[20] <- 20
#'  despike_3D(y)
#' }
#'
#' @export
despike_3D <- function(Yt, c1=2.5, c2=4, basis=NULL){
  qreg <- despike_3D.qreg(Yt, basis)
  despike_3D.interpolate(qreg, c1, c2)
}