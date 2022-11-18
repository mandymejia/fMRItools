#' despike_3D from AFNI step 1: residuals
#' 
#' Compute the quantile regression that \code{\link{despike_3D}} is based on.
#' 
#' @param Yt The data
#' 
#' @importFrom fda getbasismatrix create.fourier.basis
#' @importFrom quantreg rq
#' @importFrom stats resid
#' 
#' @return the quantile regression
#' 
#' @keywords internal
despike_3D.qreg <- function(Yt){
  nT <- length(Yt)
  basis <- getbasismatrix(
    seq(1, nT),
    create.fourier.basis(
      rangeval = c(0,nT),
      nbasis = 2*round((nT/30/2))+1
    )
  )
  quantile.reg <- rq(Yt~basis-1)
}

#' despike_3D from AFNI step 2: interpolate
#' 
#' Identify and interpolate the outliers for \code{\link{despike_3D}}.
#' 
#' @param qreg the quantile regression from \code{\link{despike_3D.qreg}}
#' @param c1 spike threshold. Default: \code{2.5}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. 
#' 
#' @importFrom stats mad fitted
#' @keywords internal
despike_3D.interpolate <- function(qreg, c1=2.5, c2=4){
  qreg_resid <- resid(qreg)
  s <- qreg_resid/mad(qreg_resid)
  s <- ifelse(s > c1, 
    c1 + (c2-c1)*tanh((s-c1)/(c2-c1)),
    s
  )
  fitted(qreg) + s*mad(qreg_resid)
}

#' despike_3D from AFNI
#' 
#' Identify and interpolate outliers. 
#' 
#' @param Yt The data
#' @param c1 spike threshold. Default: \code{2.5}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. 
#' 
#' @export
despike_3D <- function(Yt, c1=2.5, c2=4){
  qreg <- despike_3D.qreg(Yt)
  despike_3D.interpolate(qreg, c1, c2)
}

#' Temporary
#' 
#' Temporary
#' 
#' @param Yt,c1,c2 temporary
#' @importFrom fda getbasismatrix create.fourier.basis
#' @importFrom quantreg rq
#' @importFrom stats resid mad fitted
#' @keywords internal
despike_3D.original <- function(Yt, c1=2.5, c2=4){
  TIME = length(Yt)
  basis = getbasismatrix(seq(1,TIME),
                          create.fourier.basis(rangeval = c(0,TIME),nbasis=2*round((TIME/30/2))+1))
  quantile.reg = rq(Yt~basis-1)
  residuals = resid(quantile.reg)
  s = residuals/mad(residuals)
  s.prime = s
  for(id in which(s > c1)){
    s.prime[id] = c1 + (c2-c1)*tanh((s[id]-c1)/(c2-c1))
  }
  Yt.new = fitted(quantile.reg) + s.prime*mad(residuals)
  return(Yt.new)
}

