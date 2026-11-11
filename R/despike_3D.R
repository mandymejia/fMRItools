#' 3dDespike from AFNI, step 0: reference functions
#' 
#' The smooth curve that \code{\link{despike_3D}} fits to a time series, as in
#'  AFNI's \code{3dDespike}: a constant, a linear term, a quadratic term, and
#'  \code{corder} pairs of sines and cosines with period \code{nT}.
#' 
#' @param nT The number of timepoints.
#' @param corder The number of sine/cosine pairs. Default: \code{round(nT/30)},
#'  as in AFNI.
#' 
#' @return The \eqn{nT \times (3 + 2 \times corder)} matrix of reference
#'  functions.
#' 
#' @keywords internal
despike_3D.basis <- function(nT, corder=round(nT/30)){
  stopifnot(is.numeric(nT), length(nT)==1, nT >= 1)
  stopifnot(is.numeric(corder), length(corder)==1, corder >= 0)

  # Constant, linear, and quadratic terms (AFNI scales time to [-1, 1]).
  tm <- 0.5 * (nT - 1)
  tt <- (seq_len(nT) - 1 - tm) * (2 / nT)
  B <- cbind(1, tt, tt^2)

  # Sine and cosine terms, k = 1, ..., corder, at t = 0, ..., nT-1.
  if (corder > 0) {
    ang <- 2 * pi * outer(seq_len(nT) - 1, seq_len(corder)) / nT
    B <- cbind(B, sin(ang), cos(ang))
  }

  unname(B)
}

#' 3dDespike from AFNI, step 1 (L1 fit; AFNI's "OLD" method)
#' 
#' Fit the reference functions to the data with L1 (median) regression.
#' 
#' @param Yt The data vector.
#' @param basis Optional precomputed \code{\link{despike_3D.basis}}. Supply
#'  this when despiking many data vectors of equal length.
#' 
#' @return the quantile regression
#' @keywords internal
despike_3D.qreg <- function(Yt, basis=NULL){

  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package \"quantreg\" needed. Please install it", call. = FALSE)
  }

  if (is.null(basis)) { basis <- despike_3D.basis(length(Yt)) }
  stopifnot(nrow(basis) == length(Yt))
  quantile.reg <- quantreg::rq(Yt~basis-1)
}

#' 3dDespike from AFNI, step 1 (least-squares fit; AFNI's "NEW" method)
#' 
#' First, outliers are removed in a simple way: each point is compared to the
#'  median of the 9 values in the window containing it (the window is kept
#'  within the series at the ends), and replaced by that median if it differs
#'  by more than 6.789 times the median over time of the windows' MADs. Then
#'  the reference functions are fit to the modified data by least squares.
#' 
#' @param Yt The data vector.
#' @param basis The \code{\link{despike_3D.basis}}.
#' 
#' @return The fitted values.
#' @keywords internal
despike_3D.fit_NEW <- function(Yt, basis){
  nT <- length(Yt)
  stopifnot(nrow(basis) == nT)

  z <- Yt
  if (nT >= 9) {
    # Median and MAD of each window of 9. Row `r` of `W` is the window that
    #  starts at time point `r`.
    W <- stats::embed(Yt, 9)
    med_w <- matrixStats::rowMedians(W)
    mad_w <- matrixStats::rowMedians(abs(W - med_w))

    # The window for time point `i` is centered on it, but kept within the series.
    start <- pmin(pmax(seq_len(nT) - 5, 0), nT - 9)
    zme <- med_w[start + 1]
    zma <- mad_w[start + 1]

    thr <- 6.789 * stats::median(zma)
    if (thr > 0) {
      spike <- abs(Yt - zme) > thr
      z[spike] <- zme[spike]
    }
  }

  as.numeric(stats::lm.fit(basis, z)$fitted.values)
}

#' 3dDespike from AFNI, step 2
#' 
#' Identify and interpolate the outliers for \code{\link{despike_3D}}.
#' 
#' The residuals' standard deviation is estimated as
#'  \code{sqrt(pi/2) * median(abs(residuals))}, as in AFNI, and
#'  \code{s = residual / sigma}. By default, values with \code{|s| > c1} are
#'  replaced by the fitted value plus the compressed deviation
#'  \code{sign(s) * (c1 + (c2-c1)*tanh((|s|-c1)/(c2-c1))) * sigma}. (The AFNI
#'  documentation describes only \code{s > c1}, but the AFNI code edits both
#'  tails.) Only edited values change. If the residuals' scale is zero, the
#'  data are returned unchanged.
#' 
#' @param Yt The data vector.
#' @param fit The fitted values.
#' @param c1 spike threshold. Default: \code{2.5}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. 
#' @param localedit Instead, replace each value with \code{|s| >= c2} with the
#'  average of the nearest values on each side with \code{|s| < c2}. (\code{c1}
#'  is not used.) Default: \code{FALSE}.
#' 
#' @return The despiked data vector.
#' 
#' @importFrom stats median
#' @keywords internal
despike_3D.interpolate <- function(Yt, fit, c1=2.5, c2=4, localedit=FALSE){
  resid <- Yt - fit
  fsig <- sqrt(pi/2) * stats::median(abs(resid))
  if (!is.finite(fsig) || fsig <= 0) { return(Yt) }
  s <- resid / fsig

  out <- Yt
  if (!localedit) {
    c21 <- c2 - c1
    th <- function(x) { y <- tanh(x); y[x > 7] <- 1; y[x < -7] <- -1; y } # as AFNI
    hi <- s > c1
    lo <- s < -c1
    out[hi] <- fit[hi] + ( c1 + c21 * th((s[hi] - c1) / c21)) * fsig
    out[lo] <- fit[lo] + (-c1 + c21 * th((s[lo] + c1) / c21)) * fsig
  } else {
    ok <- abs(s) < c2
    idx <- seq_along(s)
    for (iv in idx[!ok]) {
      above <- idx[ok & idx > iv]
      below <- idx[ok & idx < iv]
      out[iv] <- if (length(above) > 0 && length(below) > 0) {
        0.5 * (Yt[above[1]] + Yt[below[length(below)]])
      } else if (length(above) > 0) {
        Yt[above[1]]
      } else if (length(below) > 0) {
        Yt[below[length(below)]]
      } else {
        fit[iv]
      }
    }
  }

  out
}

#' 3dDespike from AFNI
#' 
#' Identify and interpolate outliers. See 
#'  \href{https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDespike.html}{the
#'  AFNI documentation for 3dDespike} for additional information.
#' 
#' Follows \code{3dDespike}: a smooth curve is fit to the time series
#'  (\code{\link{despike_3D.basis}}), the residuals' scale is estimated, and
#'  values far from the curve are pulled in (\code{\link{despike_3D.interpolate}}).
#'  Not implemented: AFNI's \code{-NEW25} and \code{-ignore} options, and
#'  voxel masking.
#' 
#' @param Yt The data vector (at least 15 time points).
#' @param c1 spike threshold (AFNI's \code{-cut}). Default: \code{2.5}. Must be
#'  at least \code{1}.
#' @param c2 upper range of the acceptable deviation from the fit. Default: 
#'  \code{4}. Must be at least \code{c1 + 0.5}.
#' @param basis Optional precomputed \code{\link{despike_3D.basis}}, to avoid
#'  recomputing it for many data vectors of equal length.
#' @param method How to fit the curve: \code{"OLD"} for the L1 fit (requires
#'  the \code{quantreg} package), \code{"NEW"} for the faster least-squares fit
#'  after a simple outlier removal (the results are similar, not identical), or
#'  \code{"auto"} (default) to use \code{"NEW"} if there are more than 500 time
#'  points and \code{"OLD"} otherwise, as AFNI does.
#' @param localedit Use AFNI's \code{-localedit} method to replace spikes?
#'  Default: \code{FALSE}.
#' @param corder The number of sine/cosine pairs in the curve (AFNI's
#'  \code{-corder}). Default: \code{NULL} to use \code{round(length(Yt)/30)}
#'  (at most 50 for the \code{"OLD"} method). Ignored if \code{basis} is
#'  provided.
#' 
#' @return The despiked data vector.
#'
#' @references
#'  Cox RW (1996). AFNI: software for analysis and visualization of functional
#'  magnetic resonance neuroimages. \emph{Computers and Biomedical Research},
#'  29(3), 162-173. \doi{10.1006/cbmr.1996.0014}
#'
#'  Cox RW, Hyde JS (1997). Software tools for analysis and visualization of
#'  fMRI data. \emph{NMR in Biomedicine}, 10(4-5), 171-178.
#'  \doi{10.1002/(SICI)1099-1492(199706/08)10:4/5<171::AID-NBM453>3.0.CO;2-L}
#'
#' @examples
#' y <- rnorm(99) + cos(seq(99)/15)*3
#' y[20] <- 20
#' despike_3D(y, method="NEW")
#' if (requireNamespace("quantreg", quietly=TRUE)) {
#'  despike_3D(y)
#' }
#'
#' @export
despike_3D <- function(
  Yt, c1=2.5, c2=4, basis=NULL,
  method=c("auto", "OLD", "NEW"), localedit=FALSE, corder=NULL){

  stopifnot(is.numeric(Yt), !anyNA(Yt))
  Yt <- as.numeric(Yt)
  nT <- length(Yt)
  if (nT < 15) { stop("Need at least 15 time points to despike.") }
  stopifnot(is_1(c1), is_1(c2), is_1(localedit, "logical"))
  if (c1 < 1 || c2 < c1 + 0.5) { stop("Require `c1 >= 1` and `c2 >= c1 + 0.5`.") }
  method <- match.arg(method)
  if (method == "auto") { method <- if (nT > 500) { "NEW" } else { "OLD" } }

  if (is.null(basis)) {
    if (is.null(corder)) {
      corder <- round(nT / 30)
      if (corder > 50 && method == "OLD") { corder <- 50 }
    } else {
      stopifnot(is_1(corder), corder >= 0, 4*corder + 2 <= nT)
    }
    basis <- despike_3D.basis(nT, corder)
  }
  stopifnot(nrow(basis) == nT)

  fit <- if (method == "NEW") {
    despike_3D.fit_NEW(Yt, basis)
  } else {
    as.numeric(stats::fitted(despike_3D.qreg(Yt, basis)))
  }

  despike_3D.interpolate(Yt, fit, c1, c2, localedit)
}
