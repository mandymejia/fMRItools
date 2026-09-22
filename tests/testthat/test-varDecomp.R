library(testthat)

# ---- Helpers ---------------------------------------------------------------
FIELDS <- c("grand_mean", "SST", "SSB", "SSM", "SSR")

# Independent reference for ONE location; Y is M x n (complete cases only).
ref_ss <- function(Y) {
  M <- nrow(Y); n <- ncol(Y)
  if (n == 0) return(c(nS = 0, grand_mean = NaN, SST = 0, SSB = 0, SSM = 0, SSR = 0))
  gm <- mean(Y); rm <- rowMeans(Y); cm <- colMeans(Y)
  resid <- Y - outer(rm, rep(1, n)) - outer(rep(1, M), cm) + gm
  c(nS = n, grand_mean = gm,
    SST = sum((Y - gm)^2), SSB = M * sum((cm - gm)^2),
    SSM = n * sum((rm - gm)^2), SSR = sum(resid^2))
}
ref_all <- function(x) {
  M <- dim(x)[1]
  vapply(seq_len(dim(x)[3]), function(v) {
    Y <- matrix(x[, , v], nrow = M)
    ref_ss(Y[, colSums(is.na(Y)) == 0, drop = FALSE])
  }, numeric(6))
}

make_dr <- function(M = 2, N = 12, V = 6, seed = 1) {
  set.seed(seed)
  subj  <- matrix(rnorm(N * V, sd = 2), N, V)
  visit <- matrix(rnorm(M * V, sd = 0.5), M, V)
  x <- array(NA_real_, c(M, N, V))
  for (m in seq_len(M))
    x[m, , ] <- 10 + subj + rep(visit[m, ], each = N) + matrix(rnorm(N * V), N, V)
  x
}

expect_vd_equal <- function(a, b, tolerance = 1e-8) {
  expect_equal(unname(a$nS), unname(b$nS))
  for (nm in FIELDS)
    expect_equal(unname(a[[nm]]), unname(b[[nm]]), tolerance = tolerance, info = nm)
}
pick <- function(vd, i) lapply(vd[c("nS", FIELDS)], function(z) unname(z)[i])

expect_vd_matches_ref <- function(x) {
  ref <- ref_all(x)
  vd <- var_decomp(x)
  expect_equal(unname(vd$nS), unname(ref["nS", ]))
  for (nm in FIELDS) expect_equal(unname(vd[[nm]]), unname(ref[nm, ]), info = nm)
  invisible(vd)
}

# ---- var_decomp ------------------------------------------------------------
test_that("complete data matches reference for M = 2, 3, 5", {
  for (M in c(2, 3, 5)) expect_vd_matches_ref(make_dr(M = M))
})

test_that("scattered NAs: per-location complete-case results match reference", {
  x <- make_dr(N = 20, V = 10)
  set.seed(2); x[sample(length(x), 60)] <- NA
  vd <- expect_vd_matches_ref(x)
  expect_gt(length(unique(vd$nS)), 1)   # make sure the test is meaningful
})

test_that("a subject missing one visit is dropped at that location only", {
  x <- make_dr(N = 10, V = 3)
  x[2, 3, 1] <- NA
  vd <- var_decomp(x)
  expect_equal(unname(vd$nS), c(9, 10, 10))
  expect_vd_equal(pick(vd, 1),   pick(var_decomp(x[, -3, 1, drop = FALSE]), 1))
  expect_vd_equal(pick(vd, 2:3), pick(var_decomp(x[, , 2:3, drop = FALSE]), 1:2))
})

test_that("each location is independent of the others", {
  x <- make_dr(N = 12, V = 5)
  set.seed(3); x[sample(length(x), 15)] <- NA
  vd <- var_decomp(x)
  for (v in 1:5)
    expect_vd_equal(pick(vd, v), pick(var_decomp(x[, , v, drop = FALSE]), 1))
})

test_that("SST = SSB + SSM + SSR and SSR >= 0 (up to roundoff)", {
  x <- make_dr(N = 15, V = 8)
  set.seed(2); x[sample(length(x), 25)] <- NA
  vd <- var_decomp(x); ok <- vd$nS > 0
  expect_equal(unname(vd$SST[ok]), unname((vd$SSB + vd$SSM + vd$SSR)[ok]))
  expect_true(all(vd$SSR >= -1e-8))
  # Visits identical => SSR should be ~0 (cancellation stress) and SSM ~ 0
  y <- make_dr(N = 10, V = 3); y[2, , ] <- y[1, , ]
  vy <- var_decomp(y)
  expect_true(all(abs(vy$SSR) < 1e-8)); expect_true(all(abs(vy$SSM) < 1e-8))
})

test_that("invariant to subject order and visit order", {
  x <- make_dr(M = 3, N = 12, V = 5)
  set.seed(4); x[sample(length(x), 20)] <- NA
  vd <- var_decomp(x)
  expect_vd_equal(vd, var_decomp(x[, sample(12), , drop = FALSE]))
  expect_vd_equal(vd, var_decomp(x[3:1, , , drop = FALSE]))
})

test_that("shift and scale behave", {
  x <- make_dr(N = 12, V = 4)
  x[c(5, 30, 41)] <- NA
  vd <- var_decomp(x)
  sh <- var_decomp(x + 1e6)
  expect_equal(unname(sh$grand_mean), unname(vd$grand_mean + 1e6))
  for (nm in c("SST", "SSB", "SSM", "SSR"))
    expect_equal(unname(sh[[nm]]), unname(vd[[nm]]), tolerance = 1e-6, info = nm)
  sc <- var_decomp(x * 3)
  for (nm in c("SST", "SSB", "SSM", "SSR"))
    expect_equal(unname(sc[[nm]]), unname(9 * vd[[nm]]), info = nm)
})

test_that("locations with nS = 0 or 1 don't error or leak", {
  x <- make_dr(N = 6, V = 4)
  x[, , 1] <- NA          # nobody
  x[, -2, 2] <- NA        # one subject
  x[1, 1, 3] <- NA        # 5 subjects
  vd <- expect_vd_matches_ref(x)
  expect_equal(unname(vd$nS), c(0, 1, 5, 6))
  expect_vd_equal(pick(vd, 4), pick(var_decomp(x[, , 4, drop = FALSE]), 1))
})

test_that("shape edge cases", {
  x <- make_dr(N = 7, V = 1)
  expect_equal(var_decomp(matrix(x, nrow = 2)), var_decomp(x))   # 2D == V=1
  x1 <- make_dr(N = 1, V = 4)                                    # N = 1
  expect_equal(unname(var_decomp(x1)$nS), rep(1, 4))
  xn <- x; dimnames(xn) <- list(c("a", "b"), paste0("s", 1:7), "loc")
  expect_vd_equal(var_decomp(xn), var_decomp(x))                 # dimnames
  xv <- make_dr(N = 5, V = 2); xv[1, 1, 1] <- NA
  expect_output(var_decomp(xv, verbose = TRUE), "detected")
})

test_that("NaN is treated as NA", {
  x <- make_dr(N = 8, V = 3)
  xa <- x; xa[1, 2, 1] <- NA
  xn <- x; xn[1, 2, 1] <- NaN
  expect_vd_equal(var_decomp(xn), var_decomp(xa))
})

# Fails until the `is.finite` edit; drop it if you'd rather Inf not count as missing.
test_that("Inf is treated as missing", {
  x <- make_dr(N = 8, V = 3)
  xa <- x; xa[1, 2, 1] <- NA
  xi <- x; xi[1, 2, 1] <- Inf
  expect_vd_equal(var_decomp(xi), var_decomp(xa))
})

# Fails until the `stopifnot(nM >= 2)` guard.
test_that("M = 1 is rejected", {
  expect_error(var_decomp(array(rnorm(15), c(1, 5, 3))))
  expect_error(estimate_prior_from_DR(array(rnorm(15), c(1, 5, 3))))
})

# ---- mean_squares ----------------------------------------------------------
test_that("mean_squares: NA where nS < 2, matches SS / df elsewhere", {
  x <- make_dr(N = 5, V = 3); x[, -1, 1] <- NA
  vd <- var_decomp(x); ms <- mean_squares(vd)
  expect_true(is.na(ms$MSB[1]) && is.na(ms$MSR[1]))
  expect_false(anyNA(ms$MSB[2:3]))
  expect_equal(unname(ms$MSB[2:3]), unname(vd$SSB[2:3] / (vd$nS[2:3] - 1)))
  expect_equal(unname(ms$MSR[2:3]), unname(vd$SSR[2:3] / ((2 - 1) * (vd$nS[2:3] - 1))))
})