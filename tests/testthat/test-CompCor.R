# Tests for the CompCor functions: `noise_PCs` (compute noise PCs),
#  `CompCor_clean` (apply: regress PCs + DCT + nuisance + intercept),
#  `CompCor` (one-shot), `CompCor_HCP_calc` / `get_NIFTI_ROI_masks` (NIFTI
#  wrappers) and `CompCor_HCP_apply` (CIFTI wrapper, needs ciftiTools).
#
# Strategy: every function is checked against an independent reference written
#  in this file from base R (median/MAD scaling + svd + qr.resid), plus
#  mathematical properties (orthogonality of PCs to the preliminary nuisance
#  design, orthogonality of cleaned data to the regressors, Frisch-Waugh
#  equivalences, invariances). Data are small and simulated, so this runs fast.
#
# Internal functions (e.g. `noise_PCs`, `hat_Q`) are reachable because testthat
#  runs in the package namespace. They are also copied to the global
#  environment so the file can be source()d from a runner script.
for (.nm in setdiff(ls(asNamespace("fMRItools"), all.names=TRUE), ls(globalenv()))) {
  assign(.nm, get(.nm, envir=asNamespace("fMRItools")), envir=globalenv())
}
rm(.nm)

# Helpers ----------------------------------------------------------------------

# Largest difference between the columns of A and B, allowing a sign flip of
#  each column (PCs are only defined up to sign).
sign_diff <- function(A, B) {
  stopifnot(identical(dim(A), dim(B)))
  max(vapply(seq_len(ncol(A)), function(j) {
    min(max(abs(A[,j] - B[,j])), max(abs(A[,j] + B[,j])))
  }, 1))
}

# Largest fraction of any column of P (T x Q) lying in the column space of X.
#  0 means every PC is orthogonal to every column of X.
span_frac <- function(P, X) {
  qx <- qr(X); Q <- qr.Q(qx)[, seq_len(qx$rank), drop=FALSE]
  max(sqrt(colSums(crossprod(Q, P)^2)) / sqrt(colSums(P^2)))
}
# ... for the design without its constant part (columns centered).
span_frac_nonconst <- function(P, X) {
  span_frac(P, scale(X, center=TRUE, scale=FALSE))
}
# Largest |correlation| between any column of P and the vector s.
#  (Use for "is absent" checks, where the value should be ~0.)
max_cor <- function(P, s) { max(abs(stats::cor(P, s))) }
# Multiple correlation between s and the columns of P jointly (sqrt of R^2 of
#  the regression of s on P and an intercept). Use for "is present" checks:
#  PCA may return any rotation within the subspace spanned by several
#  comparably strong signals, so no single PC needs to match one signal.
span_cor <- function(P, s) {
  sqrt(1 - sum(qr.resid(qr(cbind(1, P)), s)^2) / sum((s - mean(s))^2))
}

# Reference: the noise PCs of one ROI computed directly (T x V data in).
ref_scale <- function(X, center=TRUE, scale=TRUE, tol=1e-8) {
  if (center) { X <- sweep(X, 2, apply(X, 2, median)) }
  mad <- 1.4826 * apply(abs(X), 2, median)
  keep <- mad >= tol
  X <- X[, keep, drop=FALSE]
  if (scale) { X <- sweep(X, 2, mad[keep], "/") }
  X
}
ref_noise <- function(X, center=TRUE, scale=TRUE, nPC=4, prelim=NULL) {
  if (!is.null(prelim)) { X <- qr.resid(qr(prelim), X) }
  Xs <- ref_scale(X, center, scale)
  s <- svd(Xs)
  d <- s$d^2
  if (nPC < 1) { nPC <- sum(cumsum(d) < sum(d) * nPC) + 1 }
  list(
    PCs=s$u[, seq_len(nPC), drop=FALSE], var=d[seq_len(nPC)],
    vartotal=sum(d), nPC=nPC, Xs=Xs
  )
}

# Simulated noise ROI data (T x V): a few latent time courses of decreasing
#  strength, randomly loaded onto the voxels, plus small noise.
make_roi <- function(T_, V, sds=c(8, 5, 3.5, 2.5), noise_sd=.3) {
  lat <- sapply(sds, function(s) { s * rnorm(T_) })
  Y <- lat %*% matrix(rnorm(length(sds) * V), length(sds)) +
    matrix(rnorm(T_ * V, sd=noise_sd), T_)
  attr(Y, "lat") <- lat
  Y
}

# Two known signals, orthogonal to each other and to the intercept:
#  z = slow oscillation + trend (goes in the preliminary nuisance design),
#  w = fast oscillation (does not).
make_zw <- function(T_) {
  tt <- seq_len(T_)
  z <- sin(2 * pi * tt / 37) + 1.5 * tt / T_
  z <- as.numeric(scale(z))
  w <- qr.resid(qr(cbind(1, z)), cos(2 * pi * tt / 11))
  list(z=z, w=as.numeric(scale(w)))
}
# Noise ROI data with z and w planted (random loadings) strongly enough to
#  dominate the two leading PCs.
make_planted <- function(T_, V, zw) {
  X <- make_roi(T_, V)
  X + outer(zw$z, rnorm(V, sd=40)) + outer(zw$w, rnorm(V, sd=25))
}

# A small 4D volume with three noise ROIs (cortical WM, CSF, cerebellar WM).
#  Also returns the latent time courses of the WM ROI (planted into the "data"
#  voxels too, so there is shared noise for CompCor to remove).
make_vol <- function(T_=40, dims=c(10, 24, 10), seed=1) {
  set.seed(seed)
  labels <- array(0, dim=dims)
  labels[2:9, 2:8, 2:9] <- 3001     # wm_cort
  labels[2:9, 10:16, 2:9] <- 4      # csf
  labels[2:9, 18:22, 2:9] <- 7      # wm_cblm
  nmat <- matrix(rnorm(prod(dims) * T_), prod(dims), T_)  # voxels x time
  lat_wm <- NULL
  for (lab in c(3001, 4, 7)) {
    lin <- which(labels == lab)
    Y <- make_roi(T_, length(lin), noise_sd=1)
    if (lab == 3001) { lat_wm <- attr(Y, "lat") }
    nmat[lin, ] <- t(Y)
  }
  roi_data <- array(FALSE, dim=dims); roi_data[1, 1:6, 1:6] <- TRUE   # 36 voxels
  lin <- which(roi_data)
  nmat[lin, ] <- nmat[lin, ] + t(lat_wm[, 1, drop=FALSE] %*% matrix(rnorm(length(lin), 1, .2), 1))
  list(
    nii=array(nmat, dim=c(dims, T_)), labels=labels, roi_data=roi_data,
    masks=list(wm_cort=labels == 3001, csf=labels == 4, wm_cblm=labels == 7),
    lat_wm=lat_wm, T_=T_
  )
}

# Max discrepancy between two "noise" objects (PCs up to sign, variances).
noise_diff <- function(a, b) {
  stopifnot(identical(names(a$PCs), names(b$PCs)))
  max(vapply(names(a$PCs), function(r) {
    if (is.null(a$PCs[[r]]) && is.null(b$PCs[[r]])) { return(0) }
    max(
      sign_diff(a$PCs[[r]], b$PCs[[r]]),
      max(abs(a$var[[r]] - b$var[[r]]) / abs(b$var[[r]]))
    )
  }, 1))
}
# Orient a data matrix as T x V.
orient <- function(M, T_) { if (nrow(M) == T_ && ncol(M) != T_) { M } else { t(M) } }

# new_CompCor_noise ------------------------------------------------------------
test_that("new_CompCor_noise builds the common return type", {
  x <- new_CompCor_noise(PCs=list(a=1), var=list(a=2), ROI_noise=list(a=TRUE))
  expect_s3_class(x, "CompCor_noise")
  expect_equal(names(unclass(x)), c("PCs", "var", "ROI_noise"))
  expect_equal(x$PCs, list(a=1))
  expect_equal(x$var, list(a=2))
  expect_equal(x$ROI_noise, list(a=TRUE))
})

# noise_PCs: reference match ---------------------------------------------------
test_that("noise_PCs matches a direct SVD of the centered/scaled ROI data (all center/scale combinations)", {
  set.seed(1)
  T_ <- 50
  Xn <- list(wm=make_roi(T_, 30), csf=make_roi(T_, 20))   # V differs by ROI
  for (center in c(TRUE, FALSE)) { for (scale in c(TRUE, FALSE)) {
    out <- noise_PCs(Xn, center=center, scale=scale, noise_nPC=4)
    expect_equal(names(out), c("PCs", "var", "vartotal"))
    expect_equal(names(out$PCs), c("wm", "csf"))
    for (r in names(Xn)) {
      ref <- ref_noise(Xn[[r]], center, scale, 4)
      expect_equal(dim(out$PCs[[r]]), c(T_, 4))
      expect_lt(sign_diff(out$PCs[[r]], ref$PCs), 1e-8)
      expect_equal(out$var[[r]], ref$var)
      expect_equal(out$vartotal[[r]], ref$vartotal)
    }
  }}
})

test_that("noise_PCs: PCs are orthonormal, ordered, and their variances are the squared singular values", {
  set.seed(2)
  T_ <- 50
  X <- make_roi(T_, 40)
  out <- noise_PCs(list(a=X), noise_nPC=4)
  P <- out$PCs$a
  expect_lt(max(abs(crossprod(P) - diag(4))), 1e-10)
  expect_true(all(diff(out$var$a) < 0))
  expect_lte(sum(out$var$a), out$vartotal$a + 1e-8)
  # Var of PC k = squared norm of the projection of the scaled data on it.
  Xs <- ref_scale(X)
  expect_equal(unname(colSums(crossprod(Xs, P)^2)), out$var$a)
  # Total variance = sum of squares of the scaled data.
  expect_equal(out$vartotal$a, sum(Xs^2))
})

test_that("noise_PCs: invariances (rescaling, shifting, voxel order, time order)", {
  set.seed(3)
  T_ <- 50; V <- 30
  X <- make_roi(T_, V)
  base <- noise_PCs(list(a=X), noise_nPC=4)

  # Rescaling or shifting every voxel changes nothing (center=TRUE, scale=TRUE).
  Xt <- sweep(X, 2, runif(V, .5, 4), "*") + matrix(runif(V, -3, 3), T_, V, byrow=TRUE)
  alt <- noise_PCs(list(a=Xt), noise_nPC=4)
  expect_lt(sign_diff(alt$PCs$a, base$PCs$a), 1e-8)
  expect_equal(alt$var$a, base$var$a)

  # With scale=FALSE, a global rescaling by 5 leaves the PCs alone and scales
  #  the variances by 25.
  b2 <- noise_PCs(list(a=X), scale=FALSE, noise_nPC=4)
  a2 <- noise_PCs(list(a=5 * X), scale=FALSE, noise_nPC=4)
  expect_lt(sign_diff(a2$PCs$a, b2$PCs$a), 1e-8)
  expect_equal(a2$var$a, 25 * b2$var$a)

  # Voxel order does not matter.
  perm <- sample(V)
  expect_lt(sign_diff(noise_PCs(list(a=X[, perm]), noise_nPC=4)$PCs$a, base$PCs$a), 1e-8)

  # Permuting time permutes the PCs' rows (nothing in the method is temporal).
  tp <- sample(T_)
  tperm <- noise_PCs(list(a=X[tp, ]), noise_nPC=4)
  expect_lt(sign_diff(tperm$PCs$a, base$PCs$a[tp, ]), 1e-8)
  expect_equal(tperm$var$a, base$var$a)

  # Deterministic.
  expect_identical(noise_PCs(list(a=X), noise_nPC=4), noise_PCs(list(a=X), noise_nPC=4))
})

test_that("noise_PCs finds planted signals in the leading PCs (positive control)", {
  set.seed(4)
  T_ <- 50
  zw <- make_zw(T_)
  out <- noise_PCs(list(a=make_planted(T_, 30, zw), b=make_planted(T_, 20, zw)), noise_nPC=3)
  for (r in c("a", "b")) {
    expect_gt(span_cor(out$PCs[[r]], zw$z), .95)
    expect_gt(span_cor(out$PCs[[r]], zw$w), .95)
  }
})

# noise_PCs: number of PCs -------------------------------------------------------
test_that("noise_PCs `noise_nPC`: a single number, a vector, a list, or a variance proportion", {
  set.seed(5)
  T_ <- 50
  Xn <- list(wm=make_roi(T_, 30), csf=make_roi(T_, 20))

  # Single number is recycled over ROIs.
  out <- noise_PCs(Xn, noise_nPC=2)
  expect_equal(vapply(out$PCs, ncol, 1L), c(wm=2L, csf=2L))
  expect_equal(vapply(out$var, length, 1L), c(wm=2L, csf=2L))

  # Vector / list: matched by order.
  out <- noise_PCs(Xn, noise_nPC=c(2, 3))
  expect_equal(vapply(out$PCs, ncol, 1L), c(wm=2L, csf=3L))
  out <- noise_PCs(Xn, noise_nPC=list(2, 3))
  expect_equal(vapply(out$PCs, ncol, 1L), c(wm=2L, csf=3L))
  # The leading PCs do not depend on how many are requested.
  expect_lt(sign_diff(noise_PCs(Xn, noise_nPC=c(2, 3))$PCs$csf,
    noise_PCs(Xn, noise_nPC=4)$PCs$csf[, 1:3]), 1e-8)

  # Values in (0, 1): the smallest number of PCs reaching that proportion.
  for (p in c(.5, .9)) {
    out <- noise_PCs(Xn, noise_nPC=p)
    for (r in names(Xn)) {
      s <- svd(ref_scale(Xn[[r]]))$d^2
      k <- ncol(out$PCs[[r]])
      expect_gte(sum(s[seq_len(k)]), p * sum(s))              # reaches p ...
      if (k > 1) { expect_lt(sum(s[seq_len(k - 1)]), p * sum(s)) }  # ... minimally
      expect_equal(k, sum(cumsum(s) < p * sum(s)) + 1)
      expect_equal(length(out$var[[r]]), k)
    }
  }
  # Mixed: a count for one ROI and a proportion for the other.
  out <- noise_PCs(Xn, noise_nPC=list(2, .5))
  expect_equal(ncol(out$PCs$wm), 2)
  expect_equal(ncol(out$PCs$csf), ref_noise(Xn$csf, nPC=.5)$nPC)

  # Too few values for the ROIs errors. (Extra values are currently ignored
  #  silently, not an error: `noise_PCs(Xn, noise_nPC=c(1, 2, 3))` uses 1, 2.)
  Xn3 <- c(Xn, list(cblm=make_roi(T_, 25)))
  expect_error(noise_PCs(Xn3, noise_nPC=c(1, 2)))
})

# noise_PCs: constant voxels / empty ROIs -------------------------------------
test_that("noise_PCs removes constant voxels (with a warning) and errors if all are constant", {
  set.seed(6)
  T_ <- 50
  X <- make_roi(T_, 30)
  Xc <- X; Xc[, c(3, 9)] <- 7
  expect_warning(out <- noise_PCs(list(a=Xc), noise_nPC=4), "removing 2 constant")
  ref <- ref_noise(Xc, nPC=4)
  expect_lt(sign_diff(out$PCs$a, ref$PCs), 1e-8)
  # ... same PCs as the data without those voxels.
  expect_lt(sign_diff(out$PCs$a, noise_PCs(list(a=X[, -c(3, 9)]), noise_nPC=4)$PCs$a), 1e-8)

  # Without centering, a constant (non-zero) voxel is not zero-variance.
  expect_silent(noise_PCs(list(a=Xc), center=FALSE, noise_nPC=4))

  expect_error(noise_PCs(list(a=matrix(5, T_, 4)), noise_nPC=2), "zero-variance")
})

test_that("noise_PCs skips empty ROIs and keeps the other ROIs' results", {
  set.seed(7)
  T_ <- 50
  X <- make_roi(T_, 30)
  out <- noise_PCs(list(a=X, b=matrix(0, T_, 0)), noise_nPC=3)
  expect_equal(names(out$PCs), c("a", "b"))
  expect_null(out$PCs$b)
  expect_null(out$var$b)
  expect_lt(sign_diff(out$PCs$a, noise_PCs(list(a=X), noise_nPC=3)$PCs$a), 1e-12)
  # Unnamed input is fine.
  un <- noise_PCs(list(X, X), noise_nPC=3)
  expect_null(names(un$PCs))
  expect_equal(length(un$PCs), 2)
})

# noise_PCs: prelim_nuisance ------------------------------------------------------
test_that("noise_PCs `prelim_nuisance`: PCs are orthogonal to the regressed signal, and only to it", {
  set.seed(8)
  T_ <- 50
  zw <- make_zw(T_)
  Xn <- list(a=make_planted(T_, 30, zw), b=make_planted(T_, 20, zw))
  D <- cbind(1, zw$z)                            # intercept + z (not w)

  none <- noise_PCs(Xn, noise_nPC=3)
  prel <- noise_PCs(Xn, noise_nPC=3, prelim_nuisance=D)
  prel_nc <- noise_PCs(Xn, noise_nPC=3, center=FALSE, prelim_nuisance=D)

  for (r in names(Xn)) {
    # Control: without prelim, both signals are in the PCs.
    expect_gt(span_cor(none$PCs[[r]], zw$z), .95)
    expect_gt(span_cor(none$PCs[[r]], zw$w), .95)
    # With prelim, z is gone ...
    expect_lt(max_cor(prel$PCs[[r]], zw$z), 1e-8)
    expect_lt(span_frac_nonconst(prel$PCs[[r]], D), 1e-8)
    # ... but w, which was not regressed out, is still there.
    expect_gt(span_cor(prel$PCs[[r]], zw$w), .95)
    # center=FALSE: orthogonal to the whole design, intercept included.
    expect_lt(span_frac(prel_nc$PCs[[r]], D), 1e-8)
    expect_gt(span_cor(prel_nc$PCs[[r]], zw$w), .95)
  }
})

test_that("noise_PCs `prelim_nuisance`: a rank-deficient design of several regressor types", {
  set.seed(9)
  T_ <- 50
  Xn <- list(a=make_roi(T_, 30), b=make_roi(T_, 20))
  rp <- sapply(1:3, function(j) cumsum(rnorm(T_, sd=.1)))
  spikes <- diag(T_)[, c(10, 30)]
  D <- cbind(1, dct_bases(T_, 3), rp, spikes, rp[, 1], 2 * rp[, 2] - rp[, 3])
  expect_lt(qr(D)$rank, ncol(D))

  out <- noise_PCs(Xn, noise_nPC=3, prelim_nuisance=D)
  nc <- noise_PCs(Xn, noise_nPC=3, center=FALSE, prelim_nuisance=D)
  none <- noise_PCs(Xn, noise_nPC=3)
  for (r in names(Xn)) {
    expect_lt(span_frac_nonconst(out$PCs[[r]], D), 1e-6)
    expect_lt(span_frac(nc$PCs[[r]], D), 1e-6)
    # Control: PCs from the unregressed data are not orthogonal to D.
    expect_gt(span_frac_nonconst(none$PCs[[r]], D), .05)
    # The spiked timepoints are zeroed out of the data, hence of the PCs
    #  (center=FALSE: no constant is added back).
    expect_lt(max(abs(nc$PCs[[r]][c(10, 30), ])), 1e-8)
  }
})

test_that("noise_PCs `prelim_nuisance` equals regressing first and then computing PCs; design form is irrelevant", {
  set.seed(10)
  T_ <- 50
  Xn <- list(a=make_roi(T_, 30), b=make_roi(T_, 20))
  D <- cbind(1, dct_bases(T_, 3), rnorm(T_))
  out <- noise_PCs(Xn, noise_nPC=3, prelim_nuisance=D)

  Xm <- lapply(Xn, function(x) qr.resid(qr(D), x))
  man <- noise_PCs(Xm, noise_nPC=3)
  expect_lt(noise_diff(out, man), 1e-6)

  # Reference implementation, ROI by ROI.
  for (r in names(Xn)) {
    ref <- ref_noise(Xn[[r]], nPC=3, prelim=D)
    expect_lt(sign_diff(out$PCs[[r]], ref$PCs), 1e-6)
    expect_equal(out$var[[r]], ref$var, tolerance=1e-6)
  }

  # Same projection => same PCs: column order, column scaling, orthonormal basis.
  alt <- list(
    D[, c(5, 3, 1, 4, 2)],
    sweep(D, 2, c(2, .1, 5, 1, 30), "*"),
    qr.Q(qr(D))
  )
  for (A in alt) {
    expect_lt(noise_diff(noise_PCs(Xn, noise_nPC=3, prelim_nuisance=A), out), 1e-6)
  }

  # A vector or a data.frame is accepted (converted to a matrix).
  expect_lt(noise_diff(
    noise_PCs(Xn, noise_nPC=3, prelim_nuisance=rep(1, T_)),
    noise_PCs(Xn, noise_nPC=3, prelim_nuisance=matrix(1, T_, 1))
  ), 1e-12)
  expect_lt(noise_diff(
    noise_PCs(Xn, noise_nPC=3, prelim_nuisance=as.data.frame(D)), out
  ), 1e-12)

  # NULL = no preliminary regression.
  expect_identical(
    noise_PCs(Xn, noise_nPC=3, prelim_nuisance=NULL),
    noise_PCs(Xn, noise_nPC=3)
  )

  # Wrong number of rows.
  expect_error(noise_PCs(Xn, prelim_nuisance=D[-1, ]), "rows")
})

# CompCor_clean ----------------------------------------------------------------
test_that("CompCor_clean regresses the intercept and PCs: matches qr.resid / lm.fit", {
  set.seed(11)
  T_ <- 50; V <- 7
  Xn <- list(a=make_roi(T_, 30), b=make_roi(T_, 20))
  nz <- noise_PCs(Xn, noise_nPC=3)
  P <- do.call(cbind, nz$PCs)
  Y <- matrix(rnorm(T_ * V), T_, V) + 5       # not centered
  Y_orig <- Y

  out <- CompCor_clean(Y, nz)
  expect_equal(dim(out), dim(Y))
  expect_equal(out, unname(lm.fit(cbind(1, P), Y)$residuals))
  expect_equal(out, qr.resid(qr(cbind(1, P)), Y))
  # Orthogonal to the intercept and to every PC; idempotent.
  expect_lt(max(abs(crossprod(cbind(1, P), out))), 1e-8)
  expect_lt(max(abs(colMeans(out))), 1e-10)
  expect_equal(CompCor_clean(out, nz), out)
  # The input is not modified.
  expect_identical(Y, Y_orig)

  # A `CompCor_noise` object, a plain list with `PCs`, and PCs with NULL
  #  entries (empty ROIs) all give the same result.
  obj <- new_CompCor_noise(nz$PCs, nz$var, list())
  expect_equal(CompCor_clean(Y, obj), out)
  expect_equal(CompCor_clean(Y, list(PCs=nz$PCs)), out)
  expect_equal(CompCor_clean(Y, list(PCs=c(nz$PCs, list(empty=NULL)))), out)
  # Just one ROI's PCs.
  expect_equal(CompCor_clean(Y, list(PCs=nz$PCs["a"])),
    qr.resid(qr(cbind(1, nz$PCs$a)), Y))
})

test_that("CompCor_clean: DCT and `nuisance` are added to the same simultaneous regression", {
  set.seed(12)
  T_ <- 50; V <- 7
  nz <- noise_PCs(list(a=make_roi(T_, 30), b=make_roi(T_, 20)), noise_nPC=3)
  P <- do.call(cbind, nz$PCs)
  Y <- matrix(rnorm(T_ * V), T_, V)
  N <- cbind(rnorm(T_), cumsum(rnorm(T_, sd=.2)), diag(T_)[, 12])  # incl. a spike

  # DCT: 0 = none (default), else that many cosine bases.
  expect_equal(CompCor_clean(Y, nz, DCT=0), CompCor_clean(Y, nz))
  for (k in c(1, 3, 6)) {
    expect_equal(
      CompCor_clean(Y, nz, DCT=k),
      qr.resid(qr(cbind(1, P, dct_bases(T_, k))), Y)
    )
  }
  # nuisance alone, and DCT + nuisance.
  expect_equal(CompCor_clean(Y, nz, nuisance=N), qr.resid(qr(cbind(1, P, N)), Y))
  full <- CompCor_clean(Y, nz, DCT=3, nuisance=N)
  D <- cbind(1, P, dct_bases(T_, 3), N)
  expect_equal(full, qr.resid(qr(D), Y))
  expect_lt(max(abs(crossprod(D, full))), 1e-8)
  expect_lt(abs(full[12, 1]) + abs(full[12, 2]), 1e-8)   # spiked row is zeroed

  # `nuisance` = NULL, 0 and FALSE all mean "none".
  base <- CompCor_clean(Y, nz)
  expect_equal(CompCor_clean(Y, nz, nuisance=NULL), base)
  expect_equal(CompCor_clean(Y, nz, nuisance=0), base)
  expect_equal(CompCor_clean(Y, nz, nuisance=FALSE), base)

  # Column order and PC scaling do not matter (same column space).
  P2 <- lapply(nz$PCs, function(p) sweep(p, 2, runif(ncol(p), .1, 10), "*"))
  expect_equal(CompCor_clean(Y, list(PCs=rev(P2)), DCT=3, nuisance=N[, 3:1]), full)
})

test_that("CompCor_clean removes what is in the design and keeps what is orthogonal to it", {
  set.seed(13)
  T_ <- 50
  nz <- noise_PCs(list(a=make_roi(T_, 30)), noise_nPC=3)
  P <- nz$PCs$a
  dct <- dct_bases(T_, 3)
  D <- cbind(1, P, dct)
  s <- as.numeric(qr.resid(qr(D), rnorm(T_)))   # a signal orthogonal to D
  junk <- P %*% rnorm(3) + dct %*% rnorm(3) + 3
  Y <- cbind(2 * s, 2 * s + junk, junk)
  out <- CompCor_clean(Y, nz, DCT=3)
  expect_equal(out[, 1], 2 * s)
  expect_equal(out[, 2], 2 * s)
  expect_lt(max(abs(out[, 3])), 1e-8)
  # Sums of squares split exactly into fitted + residual (centered data).
  Yc <- scale(Y, scale=FALSE)
  fit <- Yc - out
  expect_equal(colSums(Yc^2), colSums(out^2) + colSums(fit^2))
})

test_that("CompCor_clean: orthogonalization. When PCs are orthogonal to the nuisance, sequential == simultaneous regression", {
  set.seed(14)
  T_ <- 50; V <- 6
  Xn <- list(a=make_roi(T_, 30), b=make_roi(T_, 20))
  Y <- matrix(rnorm(T_ * V), T_, V)
  N <- cbind(as.numeric(scale(cumsum(rnorm(T_)))), dct_bases(T_, 2))   # no intercept
  D0 <- cbind(1, N)

  # PCs computed after regressing D0 (center=FALSE) are orthogonal to all of D0.
  nz <- noise_PCs(Xn, center=FALSE, noise_nPC=3, prelim_nuisance=D0)
  P <- do.call(cbind, nz$PCs)
  expect_lt(span_frac(P, D0), 1e-8)
  simul <- CompCor_clean(Y, nz, nuisance=N)
  seq_nuis_then_pc <- qr.resid(qr(P), qr.resid(qr(D0), Y))
  seq_pc_then_nuis <- qr.resid(qr(D0), qr.resid(qr(cbind(1, P)), Y))
  expect_equal(simul, seq_nuis_then_pc)
  expect_equal(simul, seq_pc_then_nuis)

  # Frisch-Waugh holds regardless: PCs regress on the nuisance residual PCs.
  Pr <- qr.resid(qr(D0), P)
  expect_equal(simul, qr.resid(qr(Pr), qr.resid(qr(D0), Y)))

  # Control: PCs computed WITHOUT the preliminary regression are not orthogonal
  #  to the nuisance, so removing the nuisance first and the PCs second does
  #  not give the simultaneous result.
  nz0 <- noise_PCs(Xn, noise_nPC=3)
  P0 <- do.call(cbind, nz0$PCs)
  expect_gt(span_frac_nonconst(P0, D0), .05)
  simul0 <- CompCor_clean(Y, nz0, nuisance=N)
  seq0 <- qr.resid(qr(cbind(1, P0)), qr.resid(qr(D0), Y))
  expect_gt(max(abs(simul0 - seq0)), 1e-3)
  # ... and the simultaneous result is still the full projection.
  expect_equal(simul0, qr.resid(qr(cbind(D0, P0)), Y))

  # With center=TRUE the PCs are orthogonal to the NON-constant part of D0
  #  (the median re-adds a small constant), so the nuisance-first/PC-second
  #  sequence differs from the simultaneous one only through the constant.
  nzc <- noise_PCs(Xn, center=TRUE, noise_nPC=3, prelim_nuisance=D0)
  Pc <- do.call(cbind, nzc$PCs)
  expect_lt(span_frac_nonconst(Pc, D0), 1e-8)
})

test_that("CompCor_clean handles orientation, a vector, and a square matrix", {
  set.seed(15)
  T_ <- 50; V <- 7
  nz <- noise_PCs(list(a=make_roi(T_, 30)), noise_nPC=3)
  N <- cbind(rnorm(T_), rnorm(T_))
  Y <- matrix(rnorm(T_ * V), T_, V)
  ref <- CompCor_clean(Y, nz, DCT=3, nuisance=N)

  # V x T in => V x T out, same values.
  expect_equal(CompCor_clean(t(Y), nz, DCT=3, nuisance=N), t(ref))
  # A vector is a single time series.
  expect_equal(CompCor_clean(Y[, 1], nz), CompCor_clean(Y, nz)[, 1, drop=FALSE])
  # Square: assumed T x V.
  Ysq <- matrix(rnorm(T_ * T_), T_, T_)
  expect_equal(CompCor_clean(Ysq, nz), qr.resid(qr(cbind(1, nz$PCs$a)), Ysq))
})

test_that("CompCor_clean validates its arguments", {
  set.seed(16)
  T_ <- 50
  nz <- noise_PCs(list(a=make_roi(T_, 30)), noise_nPC=3)
  Y <- matrix(rnorm(T_ * 4), T_, 4)
  N <- matrix(rnorm(T_ * 2), T_, 2)

  expect_error(CompCor_clean(Y, "a"))                            # noise not a list
  expect_error(CompCor_clean(Y, list()))                         # no `PCs`
  expect_error(CompCor_clean(Y, list(PCs=list())), "any PCs")    # empty PCs
  expect_error(CompCor_clean(Y, list(PCs=list(a=NULL))), "any PCs")
  expect_error(CompCor_clean(Y[-1, ], nz), "timepoints")         # T mismatch
  expect_error(CompCor_clean(Y, nz, nuisance=N[-1, ]), "same number")
  expect_error(CompCor_clean(Y, nz, DCT=-1))
  expect_error(CompCor_clean(Y, nz, DCT=c(1, 2)))
  expect_error(CompCor_clean(Y, nz, DCT="a"))
})

# CompCor (one-shot) -----------------------------------------------------------
test_that("CompCor requires `ROI_noise` and returns the documented structure", {
  v <- make_vol()
  expect_error(CompCor(v$nii, ROI_data=NULL), "ROI_noise")

  # No data ROI: only the noise PCs.
  out <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=3)
  expect_s3_class(out, "CompCor")
  expect_null(out$data)
  expect_s3_class(out$noise, "CompCor_noise")
  expect_equal(names(unclass(out$noise)), c("PCs", "var", "ROI_noise"))
  expect_equal(names(out$noise$PCs), names(v$masks))
  expect_equal(vapply(out$noise$PCs, ncol, 1L), c(wm_cort=3L, csf=3L, wm_cblm=3L))
  expect_equal(vapply(out$noise$PCs, nrow, 1L), c(wm_cort=40L, csf=40L, wm_cblm=40L))

  # With a data ROI that is not everything: `ROI_data` is returned too.
  out <- CompCor(v$nii, ROI_data=v$roi_data, ROI_noise=v$masks, noise_nPC=3)
  expect_equal(sum(out$ROI_data), 36)
  expect_false(is.null(out$data))
  expect_equal(prod(dim(out$data)), 36 * v$T_)
})

test_that("CompCor noise PCs match the reference for every combination of center, scale, prelim_nuisance and noise_nPC", {
  v <- make_vol()
  set.seed(20)
  ero <- c(1, 1, 1)
  fmt <- format_data(
    X=v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=3, noise_erosion=ero
  )
  P <- cbind(1, dct_bases(v$T_, 2), rnorm(v$T_))
  grid <- expand.grid(
    center=c(TRUE, FALSE), scale=c(TRUE, FALSE),
    prelim=c(FALSE, TRUE), frac=c(FALSE, TRUE)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    nPC <- if (g$frac) { .6 } else { 3 }
    prelim <- if (g$prelim) { P } else { NULL }
    out <- CompCor(
      v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=nPC,
      noise_erosion=ero, center=g$center, scale=g$scale, prelim_nuisance=prelim
    )$noise
    for (r in names(v$masks)) {
      ref <- ref_noise(fmt$X_noise[[r]], g$center, g$scale, nPC, prelim)
      info <- paste0(
        r, ": center=", g$center, " scale=", g$scale,
        " prelim=", g$prelim, " frac=", g$frac
      )
      expect_equal(ncol(out$PCs[[r]]), ncol(ref$PCs), info=info)
      expect_lt(sign_diff(out$PCs[[r]], ref$PCs), 1e-6)
      expect_equal(out$var[[r]], ref$var, tolerance=1e-6, info=info)
    }
  }
})

test_that("CompCor: `noise_nPC` given by name is matched to ROIs by name (as documented)", {
  # NOTE: the documentation says named values are matched by name, not order.
  #  If this fails with the values assigned in ROI order instead, either
  #  `format_data` or the documentation needs changing.
  v <- make_vol()
  out <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks,
    noise_nPC=c(csf=4, wm_cblm=3, wm_cort=2))
  expect_equal(vapply(out$noise$PCs, ncol, 1L), c(wm_cort=2L, csf=4L, wm_cblm=3L))
})

test_that("CompCor: `noise_nPC` is matched to ROIs by order when unnamed; erosion shrinks the noise ROIs", {
  v <- make_vol()
  # By order.
  out <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=c(2, 4, 3))
  expect_equal(vapply(out$noise$PCs, ncol, 1L), c(wm_cort=2L, csf=4L, wm_cblm=3L))

  # Erosion: 0/NULL leave the masks alone; more erosion = fewer voxels, nested.
  n_vox <- function(o) vapply(o$noise$ROI_noise, sum, 1)
  e0 <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_erosion=NULL)
  e00 <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_erosion=0)
  e1 <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_erosion=1)
  e2 <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_erosion=2)
  expect_equal(unname(n_vox(e0)), unname(vapply(v$masks, sum, 1)))
  expect_equal(e00$noise$ROI_noise, e0$noise$ROI_noise)
  expect_true(all(n_vox(e1) < n_vox(e0)))
  expect_true(all(n_vox(e2) < n_vox(e1)))
  expect_true(all(n_vox(e1) > 0))
  for (r in names(v$masks)) {
    expect_true(all(e1$noise$ROI_noise[[r]] <= e0$noise$ROI_noise[[r]]))
    expect_true(all(e2$noise$ROI_noise[[r]] <= e1$noise$ROI_noise[[r]]))
  }
  # Different ROIs, different erosion.
  mix <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_erosion=c(2, 0, 1))
  expect_equal(unname(n_vox(mix)), unname(c(n_vox(e2)[1], n_vox(e0)[2], n_vox(e1)[3])))
})

test_that("CompCor: the cleaned data is CompCor_clean of the data with the noise PCs", {
  v <- make_vol()
  set.seed(21)
  N <- cbind(rnorm(v$T_), as.numeric(scale(cumsum(rnorm(v$T_)))))
  Ydat <- matrix(v$nii, prod(dim(v$labels)), v$T_)[which(v$roi_data), ]   # V x T

  for (nuis in list(NULL, N)) {
    out <- CompCor(v$nii, ROI_data=v$roi_data, ROI_noise=v$masks,
      noise_nPC=3, noise_erosion=1, nuisance=nuis)
    got <- orient(out$data, v$T_)                                          # T x V
    expected <- CompCor_clean(t(Ydat), out$noise, nuisance=nuis)
    expect_equal(dim(got), dim(expected))
    expect_lt(max(abs(got - expected)), 1e-8)
    # Cleaned data: orthogonal to intercept, PCs and nuisance.
    D <- cbind(1, do.call(cbind, out$noise$PCs), nuis)
    expect_lt(max(abs(crossprod(D, got))), 1e-6)
  }
})

test_that("CompCor reduces shared noise in the data; `prelim_nuisance` affects the PCs, not the data directly", {
  v <- make_vol()
  set.seed(22)
  Ydat <- t(matrix(v$nii, prod(dim(v$labels)), v$T_)[which(v$roi_data), ])  # T x V
  lat1 <- v$lat_wm[, 1]
  before <- apply(Ydat, 2, function(y) abs(cor(y, lat1)))
  out <- CompCor(v$nii, ROI_data=v$roi_data, ROI_noise=v$masks, noise_nPC=4)
  after <- apply(orient(out$data, v$T_), 2, function(y) abs(cor(y, lat1)))
  expect_gt(median(before), .5)
  expect_lt(median(after), .2)

  # prelim_nuisance is NOT regressed from the data (only from the noise ROIs):
  #  a signal in the data that is neither in the PCs nor in `nuisance` survives.
  z <- make_zw(v$T_)$z
  nii2 <- v$nii
  lin <- which(v$roi_data)
  nm <- matrix(nii2, prod(dim(v$labels)), v$T_)
  nm[lin, ] <- nm[lin, ] + outer(rnorm(length(lin), 10, 1), z)
  nii2 <- array(nm, dim=dim(v$nii))
  pr <- CompCor(nii2, ROI_data=v$roi_data, ROI_noise=v$masks, noise_nPC=3,
    prelim_nuisance=cbind(1, z))
  expect_gt(median(abs(cor(orient(pr$data, v$T_), z))), .5)
  # ... whereas `nuisance` does remove it.
  nu <- CompCor(nii2, ROI_data=v$roi_data, ROI_noise=v$masks, noise_nPC=3,
    nuisance=cbind(z))
  expect_lt(max(abs(cor(orient(nu$data, v$T_), z))), 1e-8)
  # PCs themselves are orthogonal to z in the prelim run.
  expect_lt(max(vapply(pr$noise$PCs, function(p) max(abs(cor(p, z))), 1)), 1e-8)

  # Deterministic; wrong number of rows errors.
  a <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=3)
  b <- CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks, noise_nPC=3)
  expect_identical(a$noise$PCs, b$noise$PCs)
  expect_error(CompCor(v$nii, ROI_data=NULL, ROI_noise=v$masks,
    prelim_nuisance=cbind(1, z)[-1, ]), "rows")
})

# get_NIFTI_ROI_masks ----------------------------------------------------------
test_that("get_NIFTI_ROI_masks uses the documented FreeSurfer label sets", {
  labs <- c(0, 2999, 3000, 3001, 4035, 4036, 5001, 5002, 5003, 7, 46, 6, 45,
    4, 5, 14, 15, 24, 31, 43, 44, 63, 62, 250, 255, 256)
  lab <- array(c(labs, rep(0, 2 * 3 * 5 - length(labs))), dim=c(2, 3, 5))
  at <- function(l) which(lab == l)[1]

  m <- get_NIFTI_ROI_masks(lab, c("wm_cort", "wm_cblm", "csf"))
  expect_equal(names(m), c("wm_cort", "wm_cblm", "csf"))
  for (x in m) { expect_true(is.logical(x)); expect_equal(dim(x), dim(lab)) }
  for (l in c(3000, 3001, 4035, 5001, 5002)) { expect_true(m$wm_cort[at(l)], info=as.character(l)) }
  for (l in c(0, 2999, 4036, 5003, 7, 4)) { expect_false(m$wm_cort[at(l)], info=as.character(l)) }
  for (l in c(7, 46)) { expect_true(m$wm_cblm[at(l)], info=as.character(l)) }
  for (l in c(0, 6, 45, 3001, 4)) { expect_false(m$wm_cblm[at(l)], info=as.character(l)) }
  for (l in c(4, 5, 14, 15, 24, 31, 43, 44, 63, 250, 255)) { expect_true(m$csf[at(l)], info=as.character(l)) }
  for (l in c(0, 62, 256, 7, 3001)) { expect_false(m$csf[at(l)], info=as.character(l)) }

  # Default and NULL: WM (cortical) and CSF only.
  expect_equal(names(get_NIFTI_ROI_masks(lab)), c("wm_cort", "csf"))
  expect_equal(get_NIFTI_ROI_masks(lab, NULL), get_NIFTI_ROI_masks(lab))
  expect_equal(get_NIFTI_ROI_masks(lab, "csf"), m["csf"])
  expect_equal(names(get_NIFTI_ROI_masks(lab, c("csf", "csf"))), "csf")

  # A custom named list of labels.
  cu <- get_NIFTI_ROI_masks(lab, list(a=c(7, 46), b=4))
  expect_equal(names(cu), c("a", "b"))
  expect_equal(cu$a, m$wm_cblm)
  expect_equal(sum(cu$b), 1)
  expect_true(cu$b[at(4)])

  # Invalid input.
  expect_error(get_NIFTI_ROI_masks(lab, "brainstem"))
  expect_error(get_NIFTI_ROI_masks(lab, 3), "recognized")
  expect_error(get_NIFTI_ROI_masks(lab, list(a="x")))
  expect_error(get_NIFTI_ROI_masks(lab[, , 1]))                  # not 3D
})

# CompCor_HCP_calc -------------------------------------------------------------
test_that("CompCor_HCP_calc equals CompCor on the same ROIs, with or without `idx` and prelim_nuisance", {
  v <- make_vol()
  set.seed(30)
  rois <- c("wm_cort", "csf", "wm_cblm")
  ero <- c(2, 1, 2)
  idx <- 6:v$T_
  P <- cbind(1, dct_bases(length(idx), 3), rnorm(length(idx)))

  out <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=ero, idx=idx)
  expect_s3_class(out, "CompCor_noise")
  expect_equal(names(unclass(out)), c("PCs", "var", "ROI_noise"))
  expect_equal(names(out$PCs), rois)
  expect_equal(vapply(out$PCs, nrow, 1L), setNames(rep(length(idx), 3), rois))

  core <- CompCor(v$nii[, , , idx], ROI_data=NULL, ROI_noise=v$masks,
    noise_nPC=3, noise_erosion=ero)$noise
  expect_lt(noise_diff(out, core), 1e-8)
  expect_equal(out$ROI_noise, core$ROI_noise)

  # With prelim_nuisance (rows = length(idx)).
  outp <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=ero,
    idx=idx, prelim_nuisance=P)
  corep <- CompCor(v$nii[, , , idx], ROI_data=NULL, ROI_noise=v$masks,
    noise_nPC=3, noise_erosion=ero, prelim_nuisance=P)$noise
  expect_lt(noise_diff(outp, corep), 1e-6)
  for (r in rois) { expect_lt(span_frac_nonconst(outp$PCs[[r]], P), 1e-6) }

  # center/scale pass through.
  nc <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=ero,
    idx=idx, center=FALSE, scale=FALSE)
  corenc <- CompCor(v$nii[, , , idx], ROI_data=NULL, ROI_noise=v$masks,
    noise_nPC=3, noise_erosion=ero, center=FALSE, scale=FALSE)$noise
  expect_lt(noise_diff(nc, corenc), 1e-8)

  # idx = NULL is all timepoints.
  expect_lt(noise_diff(
    CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=ero),
    CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=ero, idx=seq(v$T_))
  ), 1e-12)

  # The default ROIs are cortical WM and CSF.
  expect_equal(names(CompCor_HCP_calc(v$nii, v$labels, noise_nPC=2)$PCs), c("wm_cort", "csf"))
  # ... and `ROI_noise` can be a list of labels.
  lst <- CompCor_HCP_calc(v$nii, v$labels, list(a=3001, b=c(4, 7)), noise_nPC=2)
  expect_equal(names(lst$PCs), c("a", "b"))
  expect_equal(sum(lst$ROI_noise$b), sum(v$masks$csf) + sum(v$masks$wm_cblm))
})

test_that("CompCor_HCP_calc `idx` selects and orders timepoints", {
  v <- make_vol()
  rois <- c("wm_cort", "csf")
  idx <- c(4, 9, 11:35)
  sub <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, idx=idx)
  man <- CompCor_HCP_calc(v$nii[, , , idx], v$labels, rois, noise_nPC=3)
  expect_lt(noise_diff(sub, man), 1e-8)

  # Reversed idx => PCs with reversed rows (same variances).
  fw <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, idx=idx)
  bw <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, idx=rev(idx))
  for (r in rois) {
    expect_lt(sign_diff(bw$PCs[[r]], fw$PCs[[r]][rev(seq_along(idx)), ]), 1e-8)
    expect_equal(bw$var[[r]], fw$var[[r]])
  }
})

test_that("CompCor_HCP_calc validates its arguments and reports progress only if verbose", {
  v <- make_vol()
  rois <- c("wm_cort", "csf")
  expect_error(CompCor_HCP_calc(v$nii[, , , 1], v$labels, rois))               # not 4D
  expect_error(CompCor_HCP_calc(v$nii, v$labels[-1, , ], rois))                # dim mismatch
  expect_error(CompCor_HCP_calc(v$nii, v$labels, "brainstem"))
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx=c(1, v$T_ + 1)))
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx=0:5))
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx="a"))
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx=list(1:20)))
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx=6:v$T_,
    prelim_nuisance=matrix(1, 5, 1)), "rows")
  expect_error(CompCor_HCP_calc(v$nii, v$labels, rois, idx=6:v$T_, cii="x"))   # removed argument
  expect_warning(CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=2, idx=1:5), "very few")

  expect_silent(CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=2))
  expect_output(CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=2, verbose=TRUE), "Reading labels")
})

test_that("CompCor_HCP_calc reads NIFTI file paths (== arrays)", {
  skip_if_not_installed("RNifti")
  v <- make_vol()
  f_nii <- tempfile(fileext=".nii.gz"); f_lab <- tempfile(fileext=".nii.gz")
  on.exit(unlink(c(f_nii, f_lab)), add=TRUE)
  RNifti::writeNifti(v$nii, f_nii); RNifti::writeNifti(v$labels, f_lab)
  rois <- c("wm_cort", "csf", "wm_cblm")
  from_file <- CompCor_HCP_calc(f_nii, f_lab, rois, noise_nPC=3, noise_erosion=1, idx=6:40)
  from_arr <- CompCor_HCP_calc(v$nii, v$labels, rois, noise_nPC=3, noise_erosion=1, idx=6:40)
  expect_lt(noise_diff(from_file, from_arr), 1e-4)     # data is stored as float
})

# CompCor_HCP_apply ------------------------------------------------------------
test_that("CompCor_HCP_apply validates `noise` before doing anything else", {
  expect_error(CompCor_HCP_apply(NULL, "a"))
  expect_error(CompCor_HCP_apply(NULL, list()))
})

test_that("CompCor_HCP_apply (CIFTI) == CompCor_clean == manual regression", {
  skip_if_not_installed("ciftiTools")
  set.seed(40)
  T_ <- 40; nV <- 24
  nz <- noise_PCs(list(a=make_roi(T_, 30), b=make_roi(T_, 20)), noise_nPC=3)
  Yc <- matrix(rnorm(nV * T_), nV, T_)                              # vertices x time
  xii <- tryCatch(
    ciftiTools::as_xifti(cortexL=Yc[1:12, ], cortexR=Yc[13:24, ]),
    error=function(e) skip(paste("could not build a xifti:", conditionMessage(e)))
  )
  N <- cbind(rnorm(T_), rnorm(T_))
  P <- do.call(cbind, nz$PCs)

  for (ndct in c(0, 3)) {
    out <- as.matrix(CompCor_HCP_apply(xii, nz, DCT=ndct))
    X <- cbind(1, P); if (ndct > 0) { X <- cbind(X, dct_bases(T_, ndct)) }
    expect_lt(max(abs(out - t(qr.resid(qr(X), t(Yc))))), 1e-8)
    expect_lt(max(abs(out - CompCor_clean(Yc, nz, DCT=ndct))), 1e-8)
  }
  # `nuisance_too` is the `nuisance` of CompCor_clean.
  out <- as.matrix(CompCor_HCP_apply(xii, nz, DCT=2, nuisance_too=N))
  expect_lt(max(abs(out - CompCor_clean(Yc, nz, DCT=2, nuisance=N))), 1e-8)

  # `idx`: select timepoints of the CIFTI that match the noise PCs.
  T_big <- 50
  Ybig <- matrix(rnorm(nV * T_big), nV, T_big)
  xbig <- ciftiTools::as_xifti(cortexL=Ybig[1:12, ], cortexR=Ybig[13:24, ])
  idx <- 6:45                                                       # 40 timepoints
  out <- as.matrix(CompCor_HCP_apply(xbig, nz, idx=idx))
  expect_lt(max(abs(out - CompCor_clean(Ybig[, idx], nz))), 1e-8)

  # Wrong number of timepoints errors.
  expect_error(CompCor_HCP_apply(xbig, nz))
})
