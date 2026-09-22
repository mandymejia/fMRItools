library(fMRItools)
library(testthat)

# NOTE: mostly written by Claude, edited by Damon

# ---- Data-generating helpers -----------------------------------------------

# nV x nT matrix, each row ~ N(mu_v, sigma_v^2) iid across time.
# mu / sigma can be scalars (recycled) or length-nV vectors (per-voxel).
make_bold <- function(nV = 200, nT = 100, mu = 1000, sigma = 10) {
  mu_vec    <- if (length(mu) == nV) mu else rep(mu, nV)
  sigma_vec <- if (length(sigma) == nV) sigma else rep(sigma, nV)
  noise <- matrix(rnorm(nV * nT), nrow = nV, ncol = nT)
  BOLD <- mu_vec + noise * sigma_vec
  attr(BOLD, "mu_vec") <- mu_vec
  attr(BOLD, "sigma_vec") <- sigma_vec
  BOLD
}

# add a known nuisance effect: BOLD_out[v, t] = BOLD[v, t] + beta[v] * x[t]
add_nuisance_effect <- function(BOLD, x, beta) BOLD + outer(beta, x)

# add a constant to specific columns (a crude "spike"/outlier)
add_spikes <- function(BOLD, idx, magnitude = 1e6) {
  BOLD[, idx] <- BOLD[, idx] + magnitude
  BOLD
}

# =============================================================================
# 1. Basic dimension check / warning
# =============================================================================

test_that("warns when there are more timepoints than voxels", {
  BOLD <- make_bold(nV = 5, nT = 50)
  expect_warning(norm_BOLD(BOLD, scale_by = "none"), "More time points than voxels")
})

# =============================================================================
# 2. drop_first
# =============================================================================

test_that("drop_first removes exactly the first N volumes and they cannot influence the output", {
  nV <- 90; nT <- 80; nd <- 8
  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 5)
  BOLD_corrupt <- add_spikes(BOLD, seq(nd), magnitude = 1e6)

  out_drop   <- norm_BOLD(BOLD_corrupt, drop_first = nd, hpf = NULL, scale_by = "none")
  out_manual <- norm_BOLD(BOLD_corrupt[, -seq(nd), drop = FALSE],
                           drop_first = 0, hpf = NULL, scale_by = "none")

  expect_equal(ncol(out_drop), nT - nd)
  expect_equal(out_drop, out_manual, tolerance = 1e-8)
})

test_that("drop_first errors if it would leave too few timepoints", {
  BOLD <- make_bold(nV = 20, nT = 10)
  expect_error(norm_BOLD(BOLD, drop_first = 9, hpf = NULL, scale_by = "none"))
})

# =============================================================================
# 3. nuisance
# =============================================================================

test_that("nuisance regression exactly removes a known regressor", {
  nV <- 400; nT <- 120
  BOLD <- make_bold(nV = nV, nT = nT, mu = 500, sigma = 8)
  x    <- sin(seq(nT) / 5)
  beta <- rnorm(nV, mean = 20, sd = 5)
  BOLD_n <- add_nuisance_effect(BOLD, x, beta)

  out <- norm_BOLD(BOLD_n, nuisance = matrix(x, ncol = 1), hpf = NULL, scale_by = "none")

  resid_cor <- apply(out, 1, function(v) cor(v, x))
  expect_true(all(abs(resid_cor) < 1e-6))
})

test_that("a nuisance matrix with the wrong number of rows errors", {
  BOLD <- make_bold(nV = 100, nT = 50)
  bad_nuisance <- matrix(rnorm(40), ncol = 1) # 40 rows, need 50
  expect_error(norm_BOLD(BOLD, nuisance = bad_nuisance, hpf = NULL, scale_by = "none"))
})

# =============================================================================
# 4. scrub
# =============================================================================

test_that("scrub removes exactly the flagged volumes and does not distort the rest", {
  nV <- 500; nT <- 100
  BOLD <- make_bold(nV = nV, nT = nT, mu = 800, sigma = 6)
  idx  <- c(10, 45, 46, 90)
  BOLD_spiked <- add_spikes(BOLD, idx, magnitude = 1e6)

  out_scrub  <- norm_BOLD(BOLD_spiked, scrub = idx, hpf = NULL, scale_by = "none")
  out_manual <- norm_BOLD(BOLD[, -idx, drop = FALSE], hpf = NULL, scale_by = "none")

  expect_equal(ncol(out_scrub), nT - length(idx))
  # a spike regressor for an observation is algebraically equivalent to
  # dropping that observation from the regression entirely
  expect_equal(out_scrub, out_manual, tolerance = 1e-6)
})

test_that("scrub indices (given in original, pre-drop_first numbering) are adjusted correctly", {
  nV <- 300; nT <- 60; nd <- 5
  BOLD <- make_bold(nV = nV, nT = nT, mu = 700, sigma = 4)
  spike_idx_original <- nd + 3 # 3rd retained volume, in ORIGINAL indexing
  BOLD_spiked <- add_spikes(BOLD, spike_idx_original, magnitude = 1e6)

  out <- norm_BOLD(BOLD_spiked, drop_first = nd, scrub = spike_idx_original,
                    hpf = NULL, scale_by = "none")
  expected <- norm_BOLD(BOLD[, -seq(nd), drop = FALSE][, -3, drop = FALSE],
                         hpf = NULL, scale_by = "none")

  expect_equal(ncol(out), nT - nd - 1)
  expect_equal(out, expected, tolerance = 1e-6)
})

test_that("a filtered spike regressor exactly cancels an outlier's lpf-smeared contribution", {
  nV <- 200; nT <- 150; TR <- 1
  t0 <- 75 # interior volume, away from edge-truncation effects

  BOLD_clean   <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 3)
  BOLD_outlier <- add_spikes(BOLD_clean, t0, magnitude = 5000)

  out_clean   <- norm_BOLD(BOLD_clean,   TR = TR, hpf = NULL, lpf = 0.02,
                           scrub = t0, scale_by = "none")
  out_outlier <- norm_BOLD(BOLD_outlier, TR = TR, hpf = NULL, lpf = 0.02,
                           scrub = t0, scale_by = "none")

  expect_equal(ncol(out_clean), nT - 1)
  expect_equal(out_outlier, out_clean, tolerance = 1e-6)
})

test_that("sanity check: without scrubbing, the same outlier visibly contaminates neighbouring volumes under lpf", {
  # Companion to the test above -- confirms the outlier's lpf-smeared spread
  # is large enough to matter, so the exact-cancellation result isn't
  # trivially true just because the spread was negligible to begin with.
  nV <- 200; nT <- 150; TR <- 1
  t0 <- 75

  BOLD_clean   <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 3)
  BOLD_outlier <- add_spikes(BOLD_clean, t0, magnitude = 5000)

  out_clean   <- norm_BOLD(BOLD_clean,   TR = TR, hpf = NULL, lpf = 0.02, scale_by = "none")
  out_outlier <- norm_BOLD(BOLD_outlier, TR = TR, hpf = NULL, lpf = 0.02, scale_by = "none")

  neighbor_diff <- abs(out_outlier[, t0 + 2] - out_clean[, t0 + 2])
  expect_true(all(neighbor_diff > 1))
})


# =============================================================================
# 5. TR + hpf / lpf
# =============================================================================

test_that("hpf removes a slow linear drift", {
  nV <- 300; nT <- 200; TR <- 2
  t <- seq(nT)
  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 3)
  BOLD_drift <- BOLD + outer(rep(5, nV), t)

  out <- norm_BOLD(BOLD_drift, TR = TR, hpf = 0.01, scale_by = "none")
  drift_cor <- apply(out, 1, function(v) cor(v, t))
  expect_true(all(abs(drift_cor) < 0.05))
})

test_that("hpf at its default is silently disabled (with a message) when TR is missing", {
  BOLD <- make_bold(nV = 80, nT = 50)
  expect_message(
    out <- norm_BOLD(BOLD, TR = NULL, hpf = 0.01, scale_by = "none"),
    "Setting `hpf=NULL`"
  )
  out_nofilter <- norm_BOLD(BOLD, TR = NULL, hpf = NULL, scale_by = "none")
  expect_equal(out, out_nofilter, tolerance = 1e-8)
})

test_that("a non-default hpf without TR errors", {
  BOLD <- make_bold(nV = 80, nT = 50)
  expect_error(norm_BOLD(BOLD, hpf = 0.02, scale_by = "none"))
})

test_that("lpf without TR errors", {
  BOLD <- make_bold(nV = 90, nT = 50)
  expect_error(norm_BOLD(BOLD, hpf = NULL, lpf = 0.1, scale_by = "none"))
})

test_that("lpf attenuates high-frequency content", {
  nV <- 200; nT <- 200; TR <- 0.8
  t <- seq(nT)
  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 2)
  hf_signal <- (-1)^t * 50 # near-Nyquist alternating wiggle
  BOLD_hf <- BOLD + outer(rep(1, nV), hf_signal)

  out <- norm_BOLD(BOLD_hf, TR = TR, hpf = NULL, lpf = 0.1, scale_by = "none")
  hf_var_before <- var(hf_signal)
  hf_var_after  <- mean(apply(out, 1, function(v) var(v - mean(v))))
  expect_true(hf_var_after < hf_var_before)
})

# =============================================================================
# 6. center_rows / center_cols
# =============================================================================

test_that("center_rows=TRUE (default) removes each voxel's mean over time", {
  BOLD <- make_bold(nV = 300, nT = 80, mu = 1000, sigma = 5)
  out <- norm_BOLD(BOLD, hpf = NULL, center_rows = TRUE, center_cols = FALSE, scale_by = "none")
  expect_true(all(abs(rowMeans(out)) < 1e-6))
})

test_that("center_rows=FALSE preserves each voxel's original mean", {
  mu_vec <- seq(500, 1500, length.out = 300)
  BOLD <- make_bold(nV = 300, nT = 80, mu = mu_vec, sigma = 5)
  out <- norm_BOLD(BOLD, hpf = NULL, center_rows = FALSE, center_cols = FALSE, scale_by = "none")
  expect_equal(rowMeans(out), mu_vec, tolerance = 1)
})

test_that("center_cols=TRUE removes each timepoint's spatial mean", {
  BOLD <- make_bold(nV = 500, nT = 40, mu = 1000, sigma = 5)
  out <- norm_BOLD(BOLD, hpf = NULL, center_rows = TRUE, center_cols = TRUE, scale_by = "none")
  expect_true(all(abs(colMeans(out)) < 1e-6))
})

# =============================================================================
# 7. scale_by (mean / sd / none / FUN) + scale_FUN
# =============================================================================

test_that("scale_by='none' leaves data at (centered) native scale", {
  BOLD <- make_bold(nV = 200, nT = 60, mu = 1000, sigma = 10)
  out_none <- norm_BOLD(BOLD, hpf = NULL, scale_by = "none")
  expect_equal(rowMeans(out_none), rep(0, 200), tolerance = 1e-6)
})

test_that("scale_by='mean' recovers percent signal change of a known task effect", {
  nV <- 200; nT <- 150
  mu_vec <- rep(1000, nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = mu_vec, sigma = 1) # low noise, isolate effect
  task <- rep(c(0, 1), each = nT / 2)
  pct_amplitude <- 1 # 1% signal change
  BOLD_task <- BOLD + outer(mu_vec * pct_amplitude / 100, task)

  out <- norm_BOLD(BOLD_task, hpf = NULL, scale_by = "mean", scale_sm_FWHM=0)
  amp <- apply(out, 1, function(v) diff(range(tapply(v, task, mean))))
  expect_equal(amp, rep(pct_amplitude, nV), tolerance = 0.1)
})

test_that("scale_by='sd' yields approximately unit variance per voxel", {
  nV <- 330; nT <- 300
  sigma_vec <- seq(2, 10, length.out = nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = sigma_vec)
  out <- norm_BOLD(BOLD, hpf = NULL, scale_by = "sd", scale_sm_FWHM=0)
  sds <- apply(out, 1, sd)
  expect_equal(sds, rep(1, nV), tolerance = 0.15)
})

test_that("scale_by='FUN' scales by a user-supplied measure", {
  BOLD <- make_bold(nV = 150, nT = 100, mu = 1000, sigma = 20)
  out <- norm_BOLD(BOLD, hpf = NULL, scale_by = "FUN", scale_sm_FWHM=0,
                    scale_FUN = function(x) apply(abs(x), 1, max))
  row_max <- apply(abs(out), 1, max)
  expect_equal(row_max, rep(1, 150), tolerance = 1e-6)
})

test_that("scale_by='mean' errors when the estimated mean is near zero", {
  BOLD <- make_bold(nV = 100, nT = 50, mu = 0, sigma = 1)
  expect_error(norm_BOLD(BOLD, hpf = NULL, scale_by = "mean", scale_sm_FWHM=0), "zero")
})

# =============================================================================
# 8. scale_sm_FWHM (global vs none; local without a xifti)
# =============================================================================

test_that("scale_sm_FWHM=Inf applies one global scale factor to every voxel", {
  nV <- 250; nT <- 150
  mu_vec <- seq(500, 2000, length.out = nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = mu_vec, sigma = 1)

  centered   <- norm_BOLD(BOLD, hpf = NULL, scale_by = "none")
  out_global <- norm_BOLD(BOLD, hpf = NULL, scale_by = "mean", scale_sm_FWHM = Inf)

  # every element should be divided by the exact same constant
  ratio_mat <- centered / out_global
  expect_true(diff(range(ratio_mat)) < 1e-6 * abs(mean(ratio_mat)))
})

test_that("finite scale_sm_FWHM without scale_sm_xifti warns and falls back to no smoothing", {
  BOLD <- make_bold(nV = 200, nT = 60)
  expect_warning(
    out <- norm_BOLD(BOLD, hpf = NULL, scale_by = "mean", scale_sm_FWHM = 4),
    "Skipping smoothing"
  )
  out_nosmooth <- norm_BOLD(BOLD, hpf = NULL, scale_by = "mean", scale_sm_FWHM = 0)
  expect_equal(out, out_nosmooth, tolerance = 1e-6)
})

# =============================================================================
# 9. give_stats
# =============================================================================

test_that("give_stats returns mu/sd matching known generative parameters", {
  nV <- 250; nT <- 200
  mu_vec <- seq(800, 1200, length.out = nV)
  sigma_vec <- seq(5, 15, length.out = nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = mu_vec, sigma = sigma_vec)

  out <- norm_BOLD(BOLD, hpf = NULL, scale_by = "none", give_stats = TRUE)
  expect_true(all(c("BOLD", "mu", "sd") %in% names(out)))
  expect_equal(out$mu, mu_vec, tolerance = 0.5)
  expect_equal(out$sd, sigma_vec, tolerance = 1)
})

# =============================================================================
# 10. Combined-argument scenarios
# =============================================================================

test_that("combined: drop_first + scrub + hpf + scale_by='mean' produces a finite, correctly-sized result", {
  nV <- 250; nT <- 220; nd <- 10; TR <- 1.5
  mu_vec <- rep(1000, nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = mu_vec, sigma = 2)

  BOLD <- add_spikes(BOLD, seq(nd), magnitude = 1e6)             # burn-in artifact
  BOLD <- BOLD + outer(rep(3, nV), seq(nT))                      # slow drift
  spike_idx <- c(nd + 20, nd + 21, 180)
  BOLD <- add_spikes(BOLD, spike_idx, magnitude = 1e6)                # motion spikes
  task <- rep(c(0, 1), length.out = nT)
  BOLD <- BOLD + outer(mu_vec * 0.02, task)                          # 2% task effect

  out <- norm_BOLD(BOLD, drop_first = nd, scrub = spike_idx, TR = TR, hpf = 0.01,
                    scale_by = "mean", scale_sm_FWHM=0)

  expect_equal(ncol(out), nT - nd - length(spike_idx))
  expect_false(anyNA(out))
  expect_true(all(is.finite(out)))
})

test_that("combined: nuisance regression still exactly removes the regressor when scrubbing is also applied", {
  nV <- 200; nT <- 150
  BOLD <- make_bold(nV = nV, nT = nT, mu = 900, sigma = 4)
  x <- cos(seq(nT) / 7)
  beta <- rnorm(nV, 10, 2)
  BOLD <- add_nuisance_effect(BOLD, x, beta)
  spike_idx <- c(30, 75, 76)
  BOLD <- add_spikes(BOLD, spike_idx, magnitude = 1e6)

  out <- norm_BOLD(BOLD, nuisance = matrix(x, ncol = 1), scrub = spike_idx,
                    hpf = NULL, scale_by = "none")
  x_kept <- x[-spike_idx]
  resid_cor <- apply(out, 1, function(v) cor(v, x_kept))
  expect_true(all(abs(resid_cor) < 1e-5))
})

test_that("combined: center_cols + scale_by='sd' + give_stats interact sensibly", {
  nV <- 300; nT <- 150
  sigma_vec <- seq(3, 9, length.out = nV)
  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = sigma_vec)

  out <- norm_BOLD(BOLD, hpf = NULL, center_rows = TRUE, center_cols = TRUE,
                    scale_by = "sd", give_stats = TRUE, scale_sm_FWHM=0)
  expect_true(all(c("BOLD", "mu", "sd") %in% names(out)))
  expect_equal(apply(out$BOLD, 1, sd), rep(1, nV), tolerance = 0.15)
  expect_true(all(abs(colMeans(out$BOLD)) < 0.5))
})

test_that("combined: hpf + lpf band-pass removes both drift and high-frequency noise", {
  nV <- 1050; nT <- 400; TR <- 1
  hpf <- 0.005; lpf <- 0.015
  t <- seq_len(nT)

  BOLD <- make_bold(nV = nV, nT = nT, mu = 1000, sigma = 2)
  BOLD <- BOLD + outer(rep(4, nV), t)             # slow drift (near-DC)
  wiggle <- (-1)^t * 30
  BOLD <- BOLD + outer(rep(1, nV), wiggle)        # fast wiggle (Nyquist)

  out <- norm_BOLD(BOLD, TR = TR, hpf = hpf, lpf = lpf, scale_by = "none")

  expect_false(anyNA(out))
  expect_true(all(is.finite(out)))
  expect_true(all(abs(apply(out, 1, function(v) cor(v, t))) < 0.1))
  wiggle_var_before <- var(wiggle)
  wiggle_var_after  <- mean(apply(out, 1, function(v) var(v - mean(v))))
  expect_true(wiggle_var_after < wiggle_var_before)
})

test_that("comparing mean/sd scaling", {

  # [TO DO]: write quantitative testthat checks.

  if (interactive()) {
    xii <- read_cifti(ciftiTools::ciftiTools.files()$cifti["dscalar"], idx=1)
    scale_vec <- 17 * c(as.matrix(xii))
    nV <- nrow(xii)
    nT <- 410
    TR <- .72
    BOLD <- make_bold(nV=nV, nT=nT, mu=scale_vec, sigma=.5)
    xii2 <- newdata_xifti(xii, BOLD)
    testthat::expect_lt(
      max(abs(scale_vec-rowMeans(as.matrix(xii2)))), 1
    )

    # plot(xii, zlim=c(.5, 2.2), title="scale_vec")

    rowVars_from_norm <- function(FWHM, scale_by) {
      matrixStats::rowVars(norm_BOLD(
        as.matrix(xii2), scale_by=scale_by, scale_sm_xifti=xii,
        scale_sm_FWHM=FWHM, hpf=0
      ))
    }

    # Data with variable mean, constant SD + mean scaling:
    #   expect to induce higher SD where the mean was smaller,
    #   and the extremeness of this effect decreases w/ smoothing
    xii4 <- newdata_xifti(xii, cbind(
      rowVars_from_norm(0, "mean"),
      rowVars_from_norm(20, "mean"),
      rowVars_from_norm(150, "mean"),
      rowVars_from_norm(Inf, "mean")
    ))
    plot(xii4, idx=seq(4), title=c("0","20","150","Inf"), together="idx")

    # Data with variable mean, constant SD + SD scaling:
    #   expect SD always similar across locations
    xii4 <- newdata_xifti(xii, cbind(
      rowVars_from_norm(0, "sd"),
      rowVars_from_norm(20, "sd"),
      rowVars_from_norm(150, "sd"),
      rowVars_from_norm(Inf, "sd")
    ))
    plot(xii4, idx=seq(4))

    BOLD <- make_bold(nV=nV, nT=nT, mu=88, sigma=scale_vec)
    xii2 <- newdata_xifti(xii, BOLD)

    # Data with constant mean, variable SD + mean scaling:
    #   expect SD patterns unchanged
    xii4 <- newdata_xifti(xii, cbind(
      rowVars_from_norm(0, "mean"),
      rowVars_from_norm(20, "mean"),
      rowVars_from_norm(150, "mean"),
      rowVars_from_norm(Inf, "mean")
    ))
    plot(xii4, idx=seq(4))

    # Data with constant mean, variable SD + SD scaling:
    #   expect normalization, decreasing fit w/ smoothing
    xii4 <- newdata_xifti(xii, cbind(
      rowVars_from_norm(0, "sd"),
      rowVars_from_norm(20, "sd"),
      rowVars_from_norm(150, "sd"),
      rowVars_from_norm(Inf, "sd")
    ))
    plot(xii4, idx=seq(4))
  }
})