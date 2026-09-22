#' Compute variance decomposition
#'
#' Calculate the various ANOVA sums of squares for repeated measures data.
#' For each variable, subjects with any missing measurement are excluded.
#'
#' @param x The data as a 3D array: measurements by subjects by variables.
#'  (Alternatively, a matrix that is measurements by subjects, if only one
#'  variable exists.)
#' @param verbose If \code{TRUE}, display progress of algorithm. Default:
#'  \code{FALSE}.
#' @export
#' @return The variance decomposition.
var_decomp <- function(x, verbose=FALSE) {

  # Get data dimensions.
  if (verbose) { cat("\tChecking data dimensions and missing values presence.\n") }
  d <- dim(x)
  if (length(d) == 2) {
    x <- array(x, dim=c(dim(x), 1)); d <- dim(x)
  }
  stopifnot(length(d) == 3)
  nM <- d[1]
  nS0 <- d[2]  # total subjects (includes those with missing data)

  stopifnot(nM >= 2)

  # Handle missing values: for each variable, drop subjects without complete
  #   data, and count the remaining subjects.
  # Built one visit at a time to avoid an M x N x V logical array.
  has_data <- matrix(is.finite(x[1,,]), nrow=nS0)
  if (nM > 1) {
    for (mm in seq(2, nM)) { has_data <- has_data & matrix(is.finite(x[mm,,]), nrow=nS0) }
  }
  nS <- colSums(has_data)  # length V
  if (any(!has_data)) {
    if (verbose) { cat("\t`NA`s or non-finite values detected. Excluding subjects w/o complete data, per variable.\n") }
    x[rep(!has_data, each=nM)] <- NA
  }
  rm(has_data)

  # Variance decomposition
  sub_mean <- colMeans(x, na.rm=TRUE)                 # N x V (NaN if no data)
  grand_mean <- colMeans(sub_mean, na.rm=TRUE)        # V

  # old <- vector("list")
  # old$visit_mean <- apply(x, c(1,3), mean, na.rm=TRUE)    # M x V
  # old$grand_mean2 <- array(rep(grand_mean, each=nM*nS0), dim=dim(x))
  # old$SST <- apply((x - old$grand_mean2)^2, 3, sum, na.rm=TRUE)
  # old$SSW <- apply((x - rep(sub_mean, each=nM))^2, 3, sum, na.rm=TRUE)

  # More efficient (written by Claude, verified by Damon)
  if (verbose) { cat("\tCalculating means and variance decomposition.\n") }
  visit_mean <- matrix(NA_real_, nrow=nM, ncol=d[3])
  SST <- SSW <- 0
  gm_mat <- matrix(rep(grand_mean, each=nS0), nrow=nS0)   # N x V
  for (mm in seq_len(nM)) {
    xm <- matrix(x[mm,,], nrow=nS0)                       # N x V
    visit_mean[mm,] <- colMeans(xm, na.rm=TRUE)
    SST <- SST + colSums((xm - gm_mat)^2, na.rm=TRUE)
    SSW <- SSW + colSums((xm - sub_mean)^2, na.rm=TRUE)
  }

  ### Sum of squares
  SSB <- nM * colSums((sub_mean - rep(grand_mean, each=nS0))^2, na.rm=TRUE)
  SSM <- nS * colSums((visit_mean - rep(grand_mean, each=nM))^2, na.rm=TRUE)
  SSR <- SSW - SSM # RESIDUAL/ERROR

  list(
    nS = nS,
    nM = nM,
    grand_mean = grand_mean,
    SST = SST,
    SSB = SSB, # delete?
    SSW = SSW, # delete?
    SSM = SSM,
    SSR = SSR
  )
}

#' Compute mean squares from variance decomposition
#' @param vd The variance decomposition
#' @export
#' @return The mean squares
mean_squares <- function(vd){
  n <- vd$nS; v <- vd$nM
  n[n < 2] <- NA  # need >= 2 subjects
  SST <- vd$SST
  SSW <- vd$SSM + vd$SSR
  SSB <- vd$SST - SSW
  SSM <- vd$SSM
  SSR <- vd$SSR

  ### Mean squares
  MST <- SST / (v*n - 1)
  MSW <- SSW / ((v-1)*(n))
  MSB <- SSB / (n-1)
  MSM <- SSM / (v-1)
  MSR <- SSR / ((v-1)*(n-1))
  # covXY <- (MSB - MSR) / 2 # TRUE for v != 2?

  list(
    MST = MST,
    MSW = MSW,
    MSB = MSB,
    MSM = MSM,
    MSR = MSR
  )
}