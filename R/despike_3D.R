#' 3dDespike from AFNI, step 1
#' 
#' Compute the quantile regression that \code{\link{despike_3D}} is based on.
#' 
#' @param Yt The data vector.
#' 
#' @return the quantile regression
#' 
#' @keywords internal
despike_3D.qreg <- function(Yt){

  if (!requireNamespace("fda", quietly = TRUE)) {
    stop("Package \"fda\" needed. Please install it", call. = FALSE)
  }
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package \"quantreg\" needed. Please install it", call. = FALSE)
  }

  nT <- length(Yt)
  basis <- fda::getbasismatrix(
    seq(nT),
    fda::create.fourier.basis(
      rangeval = c(0,nT),
      nbasis = 2*round((nT/30/2))+1
    )
  )
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
  s <- qreg_resid/mad(qreg_resid)
  s <- ifelse(s > c1, 
    c1 + (c2-c1)*tanh((s-c1)/(c2-c1)),
    s
  )
  fitted(qreg) + s*mad(qreg_resid)
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
#' 
#' @examples
#' if (requireNamespace("fda", quietly=TRUE) && requireNamespace("quantreg", quietly=TRUE)) {
#'  y <- rnorm(99) + cos(seq(99)/15)*3
#'  y[20] <- 20
#'  despike_3D(y)
#' }
#'
#' @export
despike_3D <- function(Yt, c1=2.5, c2=4){
  qreg <- despike_3D.qreg(Yt)
  despike_3D.interpolate(qreg, c1, c2)
}