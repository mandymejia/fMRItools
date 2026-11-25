
#' Gaussian HPF matrix
#' 
#' Matrix for applying a Gaussian-weighted local linear 
#'  detrending high-pass filter.
#' 
#' @param nT The number of time points in the data.
#' @param sigma The standard deviation of the Gaussian kernel in time points.
#'
#' @return A matrix of size nT x nT that can be used to apply a Gaussian high-pass filter to data.
#'
#' @examples
#' cutoff_sec <- 2000 # filter cutoff in seconds (HCP convention)
#' tr <- .72 # in seconds
#' hp_sigma <- 0.5 * cutoff_sec / tr (FSL convention)
#' hp_mat <- gaussian_hpf_mat(nT=500, sigma=hp_sigma)
#' @export
gaussian_hpf_mat <- function(nT, sigma) {
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  half <- max(1L, as.integer(sigma * 3)) # need >=2 points to fit a line
  L <- matrix(0, nT, nT)
  for (t in seq_len(nT)) {
    tt <- max(t - half, 1):min(t + half, nT); dt <- tt - t
    w  <- exp(-0.5 * dt^2 / sigma^2)
    A <- sum(w * dt); C <- sum(w * dt^2); N <- sum(w); den <- C * N - A^2
    L[t, tt] <- w * (C - A * dt) / den # weights giving the fitted line value at t
  }
  diag(nT) - L
}

#' Gaussian LPF matrix
#' 
#' Matrix for applying a Gaussian low-pass filter.
#' 
#' @param nT The number of time points in the data.
#' @param sigma The standard deviation of the Gaussian kernel in time points.
#'
#' @return A matrix of size nT x nT that can be used to apply a Gaussian low-pass filter to data.
#'
#' @examples
#' tr <- .72 # in seconds
#' lp_sigma <- 2.8 / tr # FSL FEAT's (optional, default off) low-pass sigma
#' lp_mat <- gaussian_lpf_mat(nT=500, sigma=lp_sigma)
#' @export
gaussian_lpf_mat <- function(nT, sigma) {
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  half <- as.integer(sigma * 20) + 2
  S <- matrix(0, nT, nT)
  for (t in seq_len(nT)) {
    tt <- max(t - half, 1):min(t + half, nT)
    w  <- exp(-0.5 * (tt - t)^2 / sigma^2)
    S[t, tt] <- w / sum(w) # re-norm over points actually used
  }
  S
}

#' Filter matrix multiplication
#' 
#' Multiplies a data matrix by a filter matrix, efficiently 
#'  handling orientations of the input and output data to 
#'  avoid unnecessary transpositions.
#' 
#' @param dat The data matrix to multiply.
#' @param fmat The filter matrix to multiply with. The first and 
#'  second dimensions match the number of indices in \code{dat}. 
#' @param dat_idx_by_var Is `dat` indices by variables? Default: \code{NULL}
#'  (infer from the dimensions of \code{dat}; an error if \code{dat} is square).
#' @param out_idx_by_var Should the result be indices by variables? Default:
#'  the same orientation as \code{dat}.
#'
#' @return The filtered/multiplied data matrix.
#' @export
#'
#' @examples
#' dat <- matrix(rnorm(7000), 100, 70)
#' fmat <- gaussian_hpf_mat(100, 2)
#' # Below all yield the same result.
#' dat_filt <- fmat_mult(dat, fmat)
#' dat_filt2 <- fmat_mult(dat, fmat, dat_idx_by_var=TRUE)
#' dat_filt3 <- t(fmat_mult(t(dat), fmat))
#' dat_filt4 <- fmat_mult(t(dat), fmat, out_idx_by_var=TRUE)
#' dat_filt5 <- t(fmat_mult(dat, fmat, out_idx_by_var=FALSE))
#'
fmat_mult <- function(dat, fmat, dat_idx_by_var = NULL, out_idx_by_var = NULL) {
  nT <- nrow(fmat)
  stopifnot(nT == ncol(fmat))

  if (is.null(dat_idx_by_var)) {
    if (nrow(dat) == nT && ncol(dat) == nT) {
      stop("Ambiguous: `dat` is square. Set `dat_idx_by_var`.")
    } else if (nrow(dat) == nT) { dat_idx_by_var <- TRUE
    } else if (ncol(dat) == nT) { dat_idx_by_var <- FALSE
    } else stop("Neither dat dim matches nrow(fmat).")
  }
  if (is.null(out_idx_by_var)) { out_idx_by_var <- dat_idx_by_var }

  out <- if (dat_idx_by_var) { fmat %*% dat } else { dat %*% t(fmat) }

  # return in the requested orientation
  if (dat_idx_by_var == out_idx_by_var) { out } else { t(out) }
}

#' Band-Pass temporal filter
#'
#' Apply a band-pass temporal filter to data using Gaussian high-pass and 
#' low-pass filters. This function is designed to replicate the behavior 
#'  of the high-pass filter in FSL's \code{-bptf} command from 
#'  \code{fslmaths}, and to use a similar strategy for the low-pass filter.
#'
#' @param dat A numeric matrix of time-series data.
#' @param HP_sigma Numeric. The standard deviation of the Gaussian kernel
#'  for the high-pass filter, in time points. Default: \code{NULL} (no
#'  high-pass filtering). Following FSL, \code{sigma = 0.5 * cutoff_sec / TR}.
#  The HCP uses \code{cutoff_sec = 2000} for "detrending-like" behavior.
#  Note that the high-pass filter removes the mean of each time series;
#  FSL FEAT and the HCP add the mean back afterwards.
#' @param LP_sigma Numeric. The standard deviation of the Gaussian kernel
#'  for the low-pass filter, in time points. Default: \code{NULL} (no low-pass
#'  filtering). FSL FEAT's optional low-pass filter (off by default) used a
#'  sigma of 2.8 seconds, i.e. \code{LP_sigma = 2.8 / TR}. 
#  The HCP does not low-pass filter.
#' @param dat_idx_by_var Is `dat` indices (timepoints) by variables? 
#'  Default: \code{TRUE}.
#' @param out_idx_by_var Should the result be indices (timepoints) by variables?
#'  Default: \code{dat_idx_by_var} (same orientation as \code{dat}).
#'
#' @return A filtered numeric object of the same dimensions as \code{dat},
#'  transposed if \code{out_idx_by_var} and \code{dat_idx_by_var} differ.
#'
#' @references 
#' \href{https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/FEAT(2f)UserGuide.html}{FSL FEAT User Guide}
#'
#' @export
#' 
#' @examples
#' cutoff_sec <- 2000 # filter cutoff in seconds (HCP convention)
#' tr <- .72 # in seconds
#' hp_sigma <- 0.5 * cutoff_sec / tr
#' dat <- matrix(rnorm(1200 * 70), 1200, 70)
#' dat_filt <- fsl_bptf(dat, HP_sigma=hp_sigma)
#' 
fsl_bptf <- function(
  dat, HP_sigma=NULL, LP_sigma=NULL, 
  dat_idx_by_var = TRUE, out_idx_by_var = dat_idx_by_var) {

  # Evaluate the default now, before `dat_idx_by_var` is reassigned below.
  force(out_idx_by_var)

  nT <- if (dat_idx_by_var) { nrow(dat) } else { ncol(dat) }

  do_HPF <- !is.null(HP_sigma) && HP_sigma > 0
  do_LPF <- !is.null(LP_sigma) && LP_sigma > 0

  if (do_HPF) {
    mat <- gaussian_hpf_mat(nT, HP_sigma)
    dat <- fmat_mult(dat, mat, dat_idx_by_var, out_idx_by_var)
    dat_idx_by_var <- out_idx_by_var
  }

  if (do_LPF) {
    mat <- gaussian_lpf_mat(nT, LP_sigma)
    dat <- fmat_mult(dat, mat, dat_idx_by_var, out_idx_by_var)
  }

  if (!do_HPF && !do_LPF) {
    if (dat_idx_by_var != out_idx_by_var) {
      dat <- t(dat)
    }
  }

  dat
}

# Same result as `fsl_bptf` for the HPF, but slower and less 
# memory efficient. Kept for reference.
# bptf_loclin <- function(Y, sigma) {
#   Y <- as.matrix(Y); nT <- nrow(Y); idx <- seq_len(nT)
#   fit <- matrix(0, nT, ncol(Y))
#   for (t in idx) {
#     x <- idx - t
#     w <- exp(-x^2 / (2 * sigma^2))
#     S0 <- sum(w); S1 <- sum(w * x); S2 <- sum(w * x^2)
#     Sy <- crossprod(w, Y); Sxy <- crossprod(w * x, Y)
#     fit[t, ] <- (S2 * Sy - S1 * Sxy) / (S0 * S2 - S1^2)  # intercept = fit at t
#   }
#   Y - fit
# }