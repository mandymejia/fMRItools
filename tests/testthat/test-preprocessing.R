# Tests for nuisance regression, temporal filtering, censoring helpers, PCA,
#  despiking, scaling, masking, and small utilities. Internal functions
#  (e.g. `hat_Q`, `fft_detrend`) are reachable because testthat runs in the
#  package namespace. Data are small and simulated, so these run quickly.

# Make internal (non-exported) functions available when this file is
# source()d from the runner script rather than run inside the package.
for (.nm in setdiff(ls(asNamespace("fMRItools"), all.names=TRUE), ls(globalenv()))) {
  assign(.nm, get(.nm, envir=asNamespace("fMRItools")), envir=globalenv())
}
rm(.nm)

# k-th DCT-II basis function of length T_ (matches `dct_bases`)
dct_basis_k <- function(T_, k) cos(pi * (seq_len(T_) - 0.5) * k / T_)

# hat_matrix / hat_Q / nuisance_regression -----------------------------------
test_that("hat_matrix and hat_Q give the projection onto the design's column space", {
  set.seed(1)
  T_ <- 60
  D <- cbind(1, seq_len(T_), rnorm(T_), rnorm(T_))

  H <- hat_matrix(D)
  expect_equal(dim(H), c(T_, T_))
  expect_equal(H, t(H))
  expect_equal(H %*% H, H)
  expect_equal(H, D %*% solve(crossprod(D), t(D)))
  expect_equal(H %*% D, D)
  expect_equal(sum(diag(H)), ncol(D))

  # Rank-deficient designs: redundant columns do not change the projection.
  D2 <- cbind(D, D[,3] + D[,4], D[,1])
  expect_equal(hat_matrix(D2), H)
  Q <- hat_Q(D2)
  expect_equal(ncol(Q), 4)
  expect_equal(crossprod(Q), diag(4))

  # An all-zero design has rank zero.
  expect_equal(dim(hat_Q(matrix(0, 10, 1))), c(10, 0))
  expect_equal(hat_matrix(matrix(0, 10, 1)), matrix(0, 10, 10))

  # Invalid input.
  expect_error(hat_Q(cbind(1, c(NA, seq_len(9)))), "NA")
})

test_that("nuisance_regression removes the design's column space", {
  set.seed(2)
  T_ <- 80; V <- 7
  Y <- matrix(rnorm(T_ * V), T_, V)
  D <- cbind(1, seq_len(T_), dct_bases(T_, 3))

  R <- nuisance_regression(Y, D)
  expect_equal(dim(R), dim(Y))
  expect_equal(R, unname(lm.fit(D, Y)$residuals))
  expect_lt(max(abs(crossprod(D, R))), 1e-8)
  expect_equal(nuisance_regression(R, D), R)                   # idempotent
  expect_equal(nuisance_regression(t(Y), D), t(R))             # V x T data
  expect_equal(nuisance_regression(Y, cbind(D, D[,2])), R)     # redundant column
  expect_equal(nuisance_regression(Y[,1], D), R[,1,drop=FALSE]) # vector data
  expect_error(nuisance_regression(matrix(1, 5, 6), D), "compatible")

  # Spike regressors zero out the censored rows and are equivalent to
  #  regressing on the remaining rows only.
  S <- flags2spikes(c(5, 20), T_)
  R2 <- nuisance_regression(Y, cbind(D, S))
  expect_lt(max(abs(R2[c(5, 20),])), 1e-10)
  keep <- setdiff(seq_len(T_), c(5, 20))
  expect_equal(R2[keep,], unname(lm.fit(D[keep,], Y[keep,])$residuals))
})

# DCT bases and conversion -----------------------------------------------------
test_that("dct_bases matches the DCT-II formula and supports `idx`", {
  T_ <- 40
  B <- dct_bases(T_, 10)
  expect_equal(dim(B), c(T_, 10))
  expect_equal(B[,3], dct_basis_k(T_, 3))
  expect_equal(B[,10], dct_basis_k(T_, 10))

  G <- crossprod(B)
  expect_equal(diag(G), rep(T_/2, 10))
  expect_lt(max(abs(G - diag(diag(G)))), 1e-10)

  expect_equal(dct_bases(T_, idx=c(2, 7)), B[, c(2, 7)])
  expect_equal(dct_bases(8, 8)[,8], rep(c(2/pi, -2/pi), 4)) # last index special case

  expect_error(dct_bases(10, 11))
  expect_error(dct_bases(10, 0))
})

test_that("dct_convert converts between number of bases and Hz", {
  expect_equal(dct_convert(200, 1, n=8), 8/400)
  expect_equal(dct_convert(200, 1, f=.02), 8)
  expect_equal(dct2Hz(200, 1, 8), .02)
  expect_equal(Hz2dct(200, 1, .02), 8)
  expect_error(dct_convert(200, 1))
  expect_error(dct_convert(200, 1, n=1, f=.1))
})

# temporal_filter ----------------------------------------------------------------
test_that("temporal_filter (DCT) high-pass, low-pass, and band-pass", {
  T_ <- 200; TR <- 1
  slow <- dct_basis_k(T_, 2)
  mid <- dct_basis_k(T_, 60)
  fast <- dct_basis_k(T_, 180)
  X <- cbind(5*slow + mid, mid + fast, 3 + mid, 5*slow + mid + fast)

  # hpf = .02 Hz -> 8 bases; lpf = .4 Hz -> bases 160, ..., 200.
  expect_equal(
    unname(temporal_filter(X, TR, hpf=.02)),
    unname(cbind(mid, mid + fast, mid, mid + fast))
  )
  expect_equal(
    unname(temporal_filter(X, TR, hpf=.02, lpf=.4)),
    unname(cbind(mid, mid, mid, mid))
  )
  expect_equal(
    unname(temporal_filter(X, TR, hpf=NULL, lpf=.4)),
    unname(cbind(5*slow + mid, mid, mid, 5*slow + mid))
  )
  expect_equal(detrend(X, TR, hpf=.02), temporal_filter(X, TR, hpf=.02))
  expect_lt(max(abs(colMeans(temporal_filter(X, TR)))), 1e-10)

  # Integer X returns the bases instead.
  nm <- temporal_filter(T_, TR, hpf=.02)
  expect_equal(dim(nm), c(T_, 8))
  expect_equal(unname(nm), dct_bases(T_, 8))
  expect_true(all(grepl(" Hz$", colnames(nm))))
  expect_equal(ncol(temporal_filter(T_, TR, hpf=.02, lpf=.4)), 8 + 41)
})

test_that("temporal_filter handles edge cases and invalid input", {
  set.seed(3)
  Xr <- matrix(rnorm(20 * 3), 20, 3)

  # No filtering requested: just center. Integer X: zero bases.
  expect_equal(temporal_filter(Xr, 1, hpf=NULL, lpf=NULL), colCenter(Xr))
  expect_equal(dim(temporal_filter(20, 1, hpf=NULL, lpf=NULL)), c(20, 0))

  # A short scan for which the default HPF rounds to zero bases: only centers.
  expect_equal(temporal_filter(Xr, .72), colCenter(Xr))

  expect_error(temporal_filter(Xr, .72, hpf=.008, lpf=.009), "too low")
  expect_error(temporal_filter(Xr, 1, hpf=.2, lpf=.1))   # lpf must exceed hpf
  expect_error(temporal_filter(Xr, 1, hpf=.6))           # above Nyquist
  expect_error(temporal_filter("a", 1))
})

test_that("temporal_filter (FFT) and fft_detrend remove only low frequencies", {
  set.seed(4)
  T_ <- 128; N <- 4
  x <- matrix(rnorm(T_ * 3), T_, 3)

  out <- fft_detrend(x, N)
  expect_true(is.numeric(out))
  expect_equal(dim(out), dim(x))
  Fx <- stats::mvfft(x); Fo <- stats::mvfft(out)
  rm_bins <- c(seq_len(N), seq(T_ - N + 2, T_))
  expect_lt(max(Mod(Fo[rm_bins,])), 1e-8)
  expect_equal(Fo[-rm_bins,], Fx[-rm_bins,])
  expect_equal(fft_detrend(x, 1), x - rep(colMeans(x), each=T_))

  # hpf = .03 Hz -> N = round(.03 * 128) = 4
  expect_equal(temporal_filter(x, 1, hpf=.03, method="FFT"), out)
  expect_error(temporal_filter(x, 1, hpf=.03, lpf=.2, method="FFT"), "Not implemented")
  expect_error(temporal_filter(T_, 1, hpf=.03, method="FFT"), "Not applicable")
})

# flags2spikes -----------------------------------------------------------------
test_that("flags2spikes creates one-hot spike regressors", {
  S <- flags2spikes(c(3, 7), 10)
  expect_equal(dim(S), c(10, 2))
  expect_equal(colSums(S), c(1, 1))
  expect_equal(which(S[,1] == 1), 3L)
  expect_equal(which(S[,2] == 1), 7L)

  lg <- rep(FALSE, 10); lg[c(3, 7)] <- TRUE
  expect_equal(flags2spikes(lg, 10), S)
  expect_equal(flags2spikes(c(7, 3, 3), 10), S)   # sorted and de-duplicated

  expect_equal(dim(flags2spikes(integer(0), 10)), c(10, 0))
  expect_equal(dim(flags2spikes(rep(FALSE, 10), 10)), c(10, 0))
  expect_silent(flags2spikes(integer(0), 10))

  expect_error(flags2spikes(c(1, 11), 10))        # out of range
  expect_error(flags2spikes(0, 10))
  expect_error(flags2spikes(1.5, 10))
  expect_error(flags2spikes(c(TRUE, FALSE), 10))  # wrong length for logical
})

# expand_RPs -------------------------------------------------------------------
test_that("expand_RPs expands to 4N columns (36P with 9 signals)", {
  set.seed(5)
  RP <- matrix(rnorm(10 * 3), 10, 3)
  E <- expand_RPs(RP)
  expect_equal(dim(E), c(10, 12))
  expect_equal(E[,1:3], RP)
  expect_equal(E[1, 4:6], c(0, 0, 0))
  expect_equal(E[-1, 4:6], diff(RP))
  expect_equal(E[,7:9], RP^2)
  expect_equal(E[,10:12], E[,4:6]^2)
  expect_equal(ncol(expand_RPs(matrix(rnorm(50 * 9), 50, 9))), 36)
  expect_equal(dim(expand_RPs(as.data.frame(RP))), c(10, 12))
  expect_error(expand_RPs(matrix("a", 3, 3)))
})

# PCA --------------------------------------------------------------------------
test_that("PCA matches prcomp and returns unit-norm directions", {
  set.seed(6)
  nV <- 300; nT <- 12; Q <- 4
  X <- matrix(rnorm(nV * nT), nV) %*% diag(seq(3, 1, length.out=nT))
  p <- prcomp(X)

  out <- PCA(X, Q=Q)
  expect_equal(length(out$d), nT)
  expect_equal(out$d, p$sdev^2)
  expect_equal(dim(out$u), c(nT, Q))
  expect_equal(unname(abs(diag(crossprod(out$u, p$rotation[, seq_len(Q)])))), rep(1, Q))
  expect_equal(crossprod(out$u), diag(Q))
  expect_null(out$v)

  out3 <- PCA(X, Q=Q, Vdim=3)
  expect_equal(dim(out3$v), c(nV, 3))
  expect_equal(sqrt(colSums(out3$v^2)), rep(1, 3))
  expect_equal(crossprod(out3$v), diag(3))

  out1 <- PCA(X, Q=Q, Vdim=1)           # single direction
  expect_equal(dim(out1$v), c(nV, 1))
  expect_equal(sqrt(sum(out1$v^2)), 1)

  expect_equal(ncol(PCA(X, Q=3, Vdim="Q")$v), 3)

  # Vdim larger than Q is capped, with a warning.
  expect_warning(PCA(X, Q=2, Vdim=5), "Vdim")
  expect_equal(ncol(suppressWarnings(PCA(X, Q=2, Vdim=5))$v), 2)

  # Centering.
  expect_equal(PCA(scale(X, scale=FALSE), center=FALSE, Q=Q)$d, out$d)
  expect_error(PCA(X + 5, center=FALSE, Q=Q), "centered")

  # Other checks.
  expect_warning(PCA(X, Q=nT + 3), "Q > number")
  expect_warning(PCA(matrix(rnorm(5 * 20), 5, 20), Q=2), "time points")
  expect_error(PCA(as.data.frame(X), Q=Q))
  expect_error(PCA(X, Q=Q, Vdim="foo"))
})

test_that("PCA estimates Q with PESEL when Q is not given", {
  skip_if_not_installed("pesel")
  set.seed(7)
  X <- matrix(rnorm(300 * 3), 300, 3) %*% diag(c(10, 8, 6)) %*% matrix(rnorm(3 * 15), 3, 15) +
    matrix(rnorm(300 * 15, sd=.1), 300, 15)
  out <- suppressWarnings(PCA(X, Q_max=8))
  expect_true(ncol(out$u) %in% 2:5)
})

# despike_3D ---------------------------------------------------------------------
test_that("despike_3D.basis are as expected", {
  nT <- 300
  tm <- (nT - 1) / 2
  tt <- (seq_len(nT) - 1 - tm) * 2 / nT
  B <- despike_3D.basis(nT)                       # corder = round(300/30) = 10
  expect_equal(dim(B), c(nT, 3 + 2 * 10))
  expect_equal(B[, 1], rep(1, nT))                # constant
  expect_equal(B[, 2], tt)                        # linear, scaled to [-1, 1]
  expect_equal(B[, 3], tt^2)                      # quadratic
  expect_equal(qr(B)$rank, ncol(B))
  has_col <- function(v) any(apply(B, 2, function(b) isTRUE(all.equal(b, v))))
  expect_true(has_col(sin(2 * pi * 4 * (seq_len(nT) - 1) / nT)))   # k = 4
  expect_true(has_col(cos(2 * pi * 4 * (seq_len(nT) - 1) / nT)))

  # Default order is round(nT/30); `corder` overrides it.
  expect_equal(ncol(despike_3D.basis(20)), 3 + 2 * 1)   # round(20/30) = 1
  expect_equal(ncol(despike_3D.basis(60)), 3 + 2 * 2)
  expect_equal(ncol(despike_3D.basis(300, corder = 0)), 3)
  expect_equal(ncol(despike_3D.basis(300, corder = 3)), 3 + 2 * 3)
})

test_that("despike_3D.interpolate: sigma = sqrt(pi/2) * median(|resid|)", {
  n <- 41
  sig <- sqrt(pi / 2)                       # sigma when median(|resid|) = 1
  Yt <- rep(c(-1, 1), length.out = n)       # |resid| = 1 almost everywhere
  big <- c(1, 10, 20, 30, 31)
  Yt[big] <- c(10, 10, -10, 10, 10) * sig   # s = +/-10
  Yt[5] <- 3 * sig                          # s = 3: above c1 = 2.5, below c2 = 4
  fit <- rep(0, n)
  s <- Yt / sig
  afni_edit <- function(s, c1, c2) {
    x <- (abs(s) - c1) / (c2 - c1)
    th <- tanh(x); th[x > 7] <- 1
    sign(s) * (c1 + (c2 - c1) * th)
  }

  # Classic edit: both tails, for several (c1, c2); unedited values are untouched.
  for (cc in list(c(2.5, 4), c(3, 5), c(2, 3))) {
    out <- despike_3D.interpolate(Yt, fit, c1 = cc[1], c2 = cc[2])
    edited <- abs(s) > cc[1]
    expect_equal(out[edited], afni_edit(s[edited], cc[1], cc[2]) * sig)
    expect_identical(out[!edited], Yt[!edited])
  }

  # Local edit: the average of the nearest non-spike values on each side.
  out <- despike_3D.interpolate(Yt, fit, c1 = 2.5, c2 = 4, localedit = TRUE)
  expect_equal(out[1], Yt[2])                                  # only a non-spike after
  expect_equal(out[10], mean(Yt[c(9, 11)]))
  expect_equal(out[20], mean(Yt[c(19, 21)]))
  expect_equal(out[c(30, 31)], rep(mean(Yt[c(29, 32)]), 2))    # adjacent spikes
  expect_identical(out[-big], Yt[-big])                        # s = 3 < c2 is untouched

  # No spread in the residuals: nothing to despike.
  expect_identical(despike_3D.interpolate(rep(1, 30), rep(1, 30)), rep(1, 30))
})

test_that("despike_3D compresses spikes of either sign and leaves the rest alone", {
  set.seed(8)
  nT <- 200
  tt <- seq_len(nT)
  # Trend within AFNI's reference functions: constant + linear + a cosine.
  trend <- 5 + 0.02 * tt + 3 * cos(2 * pi * 2 * tt / nT)
  y <- trend + runif(nT, -.2, .2)               # bounded noise: only the spikes exceed c1
  spikes <- c(50, 120)
  ys <- y; ys[spikes] <- ys[spikes] + c(20, -20)

  has_quantreg <- requireNamespace("quantreg", quietly = TRUE)
  for (method in c("NEW", if (has_quantreg) "OLD")) {
    out <- despike_3D(ys, method = method)
    expect_equal(length(out), nT)
    expect_null(names(out))
    expect_lt(max(abs(out[spikes] - trend[spikes])), 1)
    expect_equal(out[-spikes], ys[-spikes], tolerance = 1e-3)
  }

  # Up to 500 time points the default is the L1 fit.
  if (has_quantreg) { expect_equal(despike_3D(ys), despike_3D(ys, method = "OLD")) }

  # Precomputed basis, `corder`, and `localedit`.
  expect_equal(despike_3D(ys, method = "NEW", basis = despike_3D.basis(nT)), despike_3D(ys, method = "NEW"))
  expect_equal(despike_3D(ys, method = "NEW", corder = 3), despike_3D(ys, method = "NEW", basis = despike_3D.basis(nT, 3)))
  out_le <- despike_3D(ys, method = "NEW", localedit = TRUE)
  expect_equal(out_le[spikes], c(mean(ys[c(49, 51)]), mean(ys[c(119, 121)])))
  expect_equal(out_le[-spikes], ys[-spikes])
})

test_that("despike_3D uses the least-squares ('NEW') fit by default beyond 500 time points", {
  set.seed(9)
  nT <- 600
  tt <- seq_len(nT)
  trend <- 5 + 0.01 * tt + 3 * cos(2 * pi * 3 * tt / nT) + sin(2 * pi * 7 * tt / nT)
  y <- trend + runif(nT, -.2, .2)
  spikes <- c(100, 300, 301, 500)
  ys <- y; ys[spikes] <- ys[spikes] + c(20, -20, 20, 20)

  # The simple outlier removal keeps the spikes from pulling the fit.
  basis <- despike_3D.basis(nT)
  fit_clean <- despike_3D.fit_NEW(y, basis)
  fit_spiked <- despike_3D.fit_NEW(ys, basis)
  expect_lt(max(abs(fit_clean - trend)), .2)
  expect_lt(max(abs(fit_spiked - fit_clean)), .2)

  out <- despike_3D(ys)
  expect_equal(out, despike_3D(ys, method = "NEW"))
  expect_lt(max(abs(out[spikes] - trend[spikes])), 1)
  expect_equal(out[-spikes], ys[-spikes])
})

test_that("despike_3D leaves constant series alone and validates its arguments", {
  for (method in c("NEW", if (requireNamespace("quantreg", quietly = TRUE)) "OLD")) {
    expect_equal(despike_3D(rep(1, 50), method = method), rep(1, 50))
  }
  expect_error(despike_3D(rnorm(14)), "15")                  # too few time points
  expect_error(despike_3D(rnorm(50), c1 = .5))              # c1 < 1
  expect_error(despike_3D(rnorm(50), c1 = 3, c2 = 3.2))     # c2 < c1 + 0.5
  expect_error(despike_3D(c(rnorm(49), NA)))
  expect_error(despike_3D(rnorm(50), corder = 20))          # 4*corder + 2 > nT
  expect_error(despike_3D(rnorm(50), method = "bogus"))
})

# rowMedians2 / scale_med --------------------------------------------------------
test_that("rowMedians2 matches apply(median)", {
  set.seed(9)
  x <- matrix(rnorm(40), 8); x[2,3] <- NA
  expect_equal(as.numeric(rowMedians2(x, na.rm=TRUE)), apply(x, 1, median, na.rm=TRUE))
  expect_equal(as.numeric(rowMedians2(x)), apply(x, 1, median))
})

test_that("scale_med centers on the median and scales by the MAD", {
  set.seed(10)
  M <- matrix(rnorm(50 * 5, 3, 2), 50, 5)
  Z <- scale_med(M)
  expect_equal(dim(Z), dim(M))
  expect_equal(apply(Z, 2, median), rep(0, 5))
  expect_equal(1.4826 * apply(abs(Z), 2, median), rep(1, 5))

  # Rows instead of columns keeps the input orientation.
  expect_equal(scale_med(t(M), doRows=TRUE), t(Z))

  # Constant columns.
  Mc <- M; Mc[,2] <- 7
  expect_warning(scale_med(Mc), "constant")
  expect_equal(dim(suppressWarnings(scale_med(Mc))), c(50, 4))
  expect_warning(scale_med(Mc, drop_const=FALSE), "set to NA")
  Zk <- suppressWarnings(scale_med(Mc, drop_const=FALSE))
  expect_equal(dim(Zk), c(50, 5))
  expect_true(all(is.na(Zk[,2])))

  # One remaining column stays a column.
  M1 <- cbind(2, 2, rnorm(30))
  expect_equal(dim(suppressWarnings(scale_med(M1))), c(30, 1))

  # All-NA column counts as constant.
  Mn <- M; Mn[,5] <- NA
  expect_equal(ncol(suppressWarnings(scale_med(Mn))), 4)

  expect_error(scale_med(matrix(1, 10, 3)), "zero-variance")
})

# scale_timeseries / scale_design_mat / pct_sig ------------------------------------
test_that("scale_timeseries scales to percent signal, sd units, or just centers", {
  set.seed(11)
  V <- 50; T_ <- 30
  B <- matrix(rnorm(V * T_, 100, 5), V, T_)
  m <- rowMeans(B)

  Zm <- scale_timeseries(B)                       # "auto" -> mean
  expect_equal(Zm, 100 * (B - m) / m)
  expect_lt(max(abs(rowMeans(Zm))), 1e-10)

  Zs <- scale_timeseries(B, "sd")
  expect_equal(apply(Zs, 1, var), rep(1, V))
  expect_lt(max(abs(rowMeans(Zs))), 1e-10)

  expect_equal(scale_timeseries(B, "none"), B - m)

  # Missing values in the sd scaling.
  B2 <- B; B2[3, 7] <- NA
  Zs2 <- scale_timeseries(B2, "sd")
  expect_true(is.na(Zs2[3, 7]))
  expect_equal(var(Zs2[3,], na.rm=TRUE), 1)

  # Orientation: T x V input is transposed (with a warning) by default;
  #  parcellated-style data (more time points than rows) needs transpose=FALSE.
  expect_warning(scale_timeseries(t(B)), "Transposing")
  expect_equal(suppressWarnings(scale_timeseries(t(B))), Zm)
  P <- matrix(rnorm(20 * 60, 100, 5), 20, 60)
  expect_silent(Zp <- scale_timeseries(P, "mean", transpose=FALSE))
  expect_equal(Zp, 100 * (P - rowMeans(P)) / rowMeans(P))
  expect_warning(scale_timeseries(P), "Transposing")

  # Errors.
  B3 <- B; B3[5,] <- 0
  expect_error(scale_timeseries(B3, "sd"), "sds")
  B4 <- matrix(rnorm(V * T_, .05, .001), V, T_)
  expect_error(scale_timeseries(B4, "mean"), "less than 0.1")
  expect_error(scale_timeseries(B, transpose="yes"))
})

test_that("scale_design_mat divides by the largest magnitude and centers", {
  set.seed(12)
  X <- cbind(1, rnorm(20), -abs(rnorm(20)) - 5, 0, seq_len(20))
  Z <- scale_design_mat(X)
  expect_equal(colMeans(Z), rep(0, 5))
  expect_equal(Z[,1], rep(0, 20))                 # intercept becomes zero
  expect_equal(Z[,3], (X[,3] - mean(X[,3])) / max(abs(X[,3])))  # sign preserved
  expect_equal(Z[,4], rep(0, 20))                 # all-zero column: no NaN
  expect_equal(scale_design_mat(t(X), doRows=TRUE), t(Z))
  expect_error(scale_design_mat(1:5))
  expect_error(scale_design_mat(cbind(1, NA)))
})

test_that("pct_sig normalizes to percent of the center", {
  set.seed(13)
  X <- matrix(runif(40, 50, 150), 10, 4, dimnames=list(NULL, letters[1:4]))
  P <- pct_sig(X)
  expect_equal(unname(apply(P, 2, median)), rep(100, 4))
  expect_equal(P, sweep(X, 2, apply(X, 2, median), "/") * 100)
  expect_equal(pct_sig(X, by="all"), X / median(X) * 100)
  expect_equal(unname(colMeans(pct_sig(X, center=mean))), rep(100, 4))
  expect_error(pct_sig(X, center="median"))
})

# Matrix/mask utilities ----------------------------------------------------------------
test_that("unvec_mat places vectors into the mask", {
  x <- cbind(seq(3), seq(2, 4), seq(3, 5))
  mask <- matrix(c(rep(TRUE, 3), FALSE), ncol=2)

  out <- unvec_mat(x, mask, 0)
  expect_equal(length(out), 3)
  expect_equal(names(out), c("1", "2", "3"))
  expect_equal(out[[1]], matrix(c(1, 2, 3, 0), 2))
  expect_equal(out[[2]], matrix(c(2, 3, 4, 0), 2))
  expect_equal(out[[3]], matrix(c(3, 4, 5, 0), 2))

  expect_true(is.na(unvec_mat(x, mask)[[1]][2, 2]))   # default fill
  expect_equal(unvec_mat(x, 1 * mask, 0), out)         # numeric mask
  expect_equal(unvec_mat(1:3, mask, 0)[[1]], out[[1]]) # plain vector

  dn <- list(c("a", "b"), c("c", "d"))
  mask_dn <- mask; dimnames(mask_dn) <- dn
  expect_equal(dimnames(unvec_mat(x, mask_dn)[[1]]), dn)

  expect_error(unvec_mat(x[1:2,], mask))               # wrong number of rows
})

test_that("unmask_mat and colCenter work", {
  x <- matrix(seq_len(6), 2, 3)
  um <- unmask_mat(x, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(dim(um), c(4, 3))
  expect_equal(um[c(1, 3),], x)
  expect_true(all(is.na(um[c(2, 4),])))
  um2 <- unmask_mat(x, c(TRUE, FALSE, TRUE, TRUE, FALSE), mask_dim=2, fill=0)
  expect_equal(dim(um2), c(2, 5))
  expect_equal(um2[,c(1, 3, 4)], x)
  expect_equal(um2[,c(2, 5)], matrix(0, 2, 2))
  expect_error(unmask_mat(x, c(TRUE, FALSE)))

  set.seed(14)
  X <- matrix(rnorm(40, 5), 10)
  expect_equal(colMeans(colCenter(X)), rep(0, 4))
  expect_equal(colCenter(X)[,2], X[,2] - mean(X[,2]))
})

test_that("mask_BOLD flags NA, low-variance, low-mean, and non-finite-variance locations", {
  set.seed(15)
  B <- matrix(rnorm(20 * 30), 20, 30)
  B[3, 5] <- NA            # NA
  B[7,] <- 1               # zero variance
  B[9,] <- rnorm(30, -100) # low mean
  B[11, 2] <- Inf          # NaN variance

  m <- mask_BOLD(B, verbose=FALSE)
  expect_true(is.logical(m)); expect_equal(length(m), 20)
  expect_equal(which(!m), c(3L, 7L, 11L))
  expect_equal(which(!mask_BOLD(B, meanTol=-50, verbose=FALSE)), c(3L, 7L, 9L, 11L))
  expect_equal(which(!mask_BOLD(B, varTol=0, verbose=FALSE)), 3L)

  expect_silent(mask_BOLD(B, verbose=FALSE))
  expect_output(mask_BOLD(B), "NA/NaN")
  expect_output(mask_BOLD(B), "low variance")

  B2 <- B; B2[,1] <- NA
  expect_true(!any(mask_BOLD(B2, verbose=FALSE)))
  expect_equal(mask_BOLD(B, verbose=FALSE), m)
  expect_error(mask_BOLD(1:5))
})

# Argument-checking utilities ----------------------------------------------------------------
test_that("match_exactly returns FALSE on mismatch and reports per `fail_action`", {
  expect_true(match_exactly(c("a", "b"), c("a", "b")))
  expect_false(match_exactly("a", "b", "nothing"))
  expect_false(match_exactly(c("a", "b"), c("b", "a"), "nothing"))
  expect_false(match_exactly("a", c("a", "b"), "nothing"))
  expect_message(match_exactly("a", "b"), "Mismatch")
  expect_warning(match_exactly("a", "b", "warning"), "Mismatch")
  expect_error(match_exactly("a", "b", "stop"), "Mismatch")
})

test_that("match_input matches user input to expected values", {
  expected <- c("left", "right", "subcortical")
  expect_equal(match_input("le", expected), "left")
  expect_equal(match_input(c("right", "sub"), expected), c("right", "subcortical"))
  expect_error(match_input("x", expected))
  expect_warning(match_input("x", expected, "warning"))
  expect_null(match_input("x", expected, "nothing"))
})

test_that("is_1, is_posNum and is_integer reject NA and non-finite values", {
  expect_true(is_1(3))
  expect_true(is_1("a", "character"))
  expect_true(is_1(TRUE, "logical"))
  expect_false(is_1(c(1, 2)))
  expect_false(is_1("a"))
  expect_false(is_1(NA))
  expect_false(is_1(NA, "logical"))
  expect_false(is_1(NA_real_))
  expect_false(is_1(NaN))

  expect_true(is_posNum(3))
  expect_false(is_posNum(0))
  expect_true(is_posNum(0, zero_ok=TRUE))
  expect_false(is_posNum(-1))
  expect_false(is_posNum(NA_real_))

  expect_true(is_integer(3))
  expect_false(is_integer(2.5))
  expect_true(is_integer(-1))
  expect_false(is_integer(-1, nneg=TRUE))
  expect_false(is_integer(NA_real_))
  expect_false(is_integer(Inf))
  expect_false(is_integer("a"))

  expect_true(all_integers(c(1, 2, 3)))
  expect_false(all_integers(c(1, 2.5)))
  expect_true(all_binary(c(0, 1, 1)))
  expect_false(all_binary(c(0, 2)))
  expect_true(all_binary(c(TRUE, FALSE)))
  expect_false(all_binary(c(TRUE, FALSE), logical_ok=FALSE))
  expect_true(is_constant(rep(2, 5)))
  expect_false(is_constant(c(1, 2)))
})

test_that("get_FORMAT maps formats and errors on unknown ones", {
  expect_equal(get_FORMAT("xifti"), "CIFTI")
  expect_equal(get_FORMAT("gifti2"), "GIFTI2")
  expect_equal(get_FORMAT("nifti"), "NIFTI")
  expect_equal(get_FORMAT("RDS"), "MATRIX")
  expect_error(get_FORMAT("bogus"), "Unrecognized")
  expect_null(check_req_ifti_pkg("MATRIX"))
})

test_that("skew_pos, sign_flip and Mode work", {
  s <- c(1, 2, 10, 1, 2)
  expect_true(skew_pos(s))
  expect_false(skew_pos(-s))
  expect_equal(unname(skew_pos(cbind(s, -s))), c(TRUE, FALSE))

  M <- matrix(seq_len(6), 3, 2)
  flipped <- sign_flip(list(S=cbind(s, -s), M=M))
  expect_true(all(skew_pos(flipped$S)))
  expect_equal(unname(flipped$S[,2]), s)
  expect_equal(flipped$M[,1], M[,1])
  expect_equal(flipped$M[,2], -M[,2])
  expect_error(sign_flip(list(S=cbind(s, -s))))

  expect_equal(Mode(c(1, 2, 2, 3)), 2)
  expect_equal(Mode(c(1, 5, 2)), 2)   # all unique: median
})

# dice_coef ----------------------------------------------------------------------------------
test_that("dice_coef computes pairwise Dice coefficients", {
  x <- c(1, 1, 0, 0); y <- c(1, 0, 1, 0)
  expect_equal(dice_coef(x, y), matrix(.5, 1, 1))
  expect_equal(dice_coef(x, x), matrix(1, 1, 1))
  expect_equal(dice_coef(x, c(0, 0, 1, 1)), matrix(0, 1, 1))

  X <- cbind(x, y)
  D <- dice_coef(X)
  expect_equal(dim(D), c(2, 2))
  expect_equal(unname(diag(D)), c(1, 1))
  expect_equal(unname(D[1, 2]), .5)
  expect_equal(unname(D), unname(t(D)))
  expect_equal(unname(dice_coef(cbind(x, 0), cbind(y, 0))[2, 2]), 0)  # empty vs empty

  expect_error(dice_coef(x, c(1, 0)))
})

# bandstop_filter ----------------------------------------------------------------------------
test_that("bandstop_filter attenuates the stopband and preserves other frequencies", {
  skip_if_not_installed("gsignal")
  TR <- .72; T_ <- 400
  tt <- (seq_len(T_) - 1) * TR
  inband <- sin(2*pi*.37*tt)   # within the .31-.43 Hz stopband
  outband <- sin(2*pi*.05*tt)
  X <- cbind(inband, outband)

  out <- bandstop_filter(X, TR, .31, .43)
  expect_equal(dim(out), dim(X))
  mid <- 80:320                # ignore edge transients
  expect_lt(sd(out[mid, 1]), .2 * sd(inband[mid]))
  expect_lt(max(abs(out[mid, 2] - outband[mid])), .1)

  expect_equal(NROW(bandstop_filter(stats::ts(inband), TR, .31, .43)), T_)
  expect_equal(NROW(bandstop_filter(inband, TR, .31, .43)), T_)

  expect_error(bandstop_filter(X, TR, .31, .70))   # above Nyquist
  expect_error(bandstop_filter(X, TR, .43, .31))   # f1 >= f2
  expect_error(bandstop_filter(rnorm(5), TR, .31, .43))
})
