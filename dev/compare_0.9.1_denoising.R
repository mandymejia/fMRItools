# =============================================================================
# Old-vs-new comparison for every function changed
# using Claude, after Commit 3057439 on 9.0 (shortly after 0.8 pushed to CRAN)
#
# Usage (paths default to locations relative to this script, which lives in dev/):
#   Rscript compare_0_9_1_0_9_0.R [new_R_dir] [old_R_dir]
# defaults: new = ../R          (current package source)
#           old = 0.9.0/R       (frozen copy of the previous release's R/)
#
# Each check is one of:
#   SAME           old and new are expected to give identical results
#                  (to a tolerance), or both to error. Rows labeled
#                  "[k settings]" run k different argument settings.
#   INTENDED_DIFF  the change fixes a bug / changes behavior on purpose; the
#                  check PASSES if the documented difference is observed (and,
#                  where noted, if the difference is confined to its scope).
#   SKIP           an optional package needed for the check is not installed.
#   COVERAGE       functions whose code differs between OLD and NEW (found
#                  automatically) but that no check above calls out (GAP), and
#                  unchanged functions that call changed ones (INFO).
#
# Old functions live in env OLD (every file in the old dir, untouched).
# New functions live in env NEW (every file in the new dir).
# Functions keep their defining environment as closure, so OLD functions call
# OLD helpers and NEW functions call NEW helpers.
# =============================================================================

script_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(o)) return(dirname(normalizePath(o)))
  getwd()
})
args <- commandArgs(trailingOnly = TRUE)
new_dir <- if (length(args) >= 1) args[1] else file.path(script_dir, "..", "R")
old_dir <- if (length(args) >= 2) args[2] else file.path(script_dir, "0.9.0", "R")
stopifnot(dir.exists(old_dir), dir.exists(new_dir))

load_dir <- function(dir) {
  env <- new.env(parent = globalenv())
  for (f in list.files(dir, pattern = "\\.[Rr]$", full.names = TRUE)) {
    tryCatch(
      sys.source(f, envir = env, keep.source = FALSE),
      error = function(e) warning("Could not source ", basename(f), ": ", conditionMessage(e), call. = FALSE)
    )
  }
  env
}
OLD <- load_dir(old_dir)
NEW <- load_dir(new_dir)
# Guard against swapped directories: hat_Q is new in the current code.
if (!is.function(NEW$hat_Q) || is.function(OLD$hat_Q)) {
  stop("`new_dir` should contain the current code (has hat_Q) and `old_dir` the 0.9.0 code (no hat_Q). ",
       "Usage: Rscript <script> <new_R_dir> <old_R_dir>. Got new = ", new_dir, ", old = ", old_dir, call. = FALSE)
}

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# What differs between OLD and NEW (formatting/comments are ignored).
fn_names <- function(env) Filter(function(n) is.function(get(n, envir = env)), ls(env, all.names = TRUE))
fn_text  <- function(env, n) paste(deparse(get(n, envir = env)), collapse = "\n")
fn_old <- fn_names(OLD); fn_new <- fn_names(NEW)
fn_added    <- setdiff(fn_new, fn_old)
fn_removed  <- setdiff(fn_old, fn_new)
fn_modified <- Filter(function(n) fn_text(OLD, n) != fn_text(NEW, n), intersect(fn_old, fn_new))

# ---- harness ----------------------------------------------------------------
results <- list()
record <- function(label, kind, status, detail = "") {
  results[[length(results) + 1]] <<- data.frame(
    label = label, kind = kind, status = status, detail = detail,
    stringsAsFactors = FALSE
  )
}
# Functions that the checks below exercise (used for the coverage report).
covered <- character()
covers <- function(...) covered <<- union(covered, c(...))

# Evaluate an expression, capturing value / error / warnings.
run <- function(expr) {
  w <- character()
  v <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
    ),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "run_error")
  )
  list(value = v, error = inherits(v, "run_error"), warnings = w)
}
errmsg <- function(r) if (r$error) r$value$msg else "none"
same <- function(a, b, tol = 1e-8) isTRUE(all.equal(a, b, tolerance = tol))
msgs_of <- function(expr) {   # text of any message() calls
  m <- character()
  withCallingHandlers(
    tryCatch(expr, error = function(e) NULL),
    message = function(x) { m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage") }
  )
  m
}

# Expect old and new to agree (or both to error).
cmp <- function(label, old_expr, new_expr, tol = 1e-8) {
  ro <- run(old_expr); rn <- run(new_expr)
  if (ro$error && rn$error) {
    record(label, "SAME", "PASS", "both error")
  } else if (ro$error || rn$error) {
    record(label, "SAME", "FAIL",
      paste0("old error: ", errmsg(ro), " | new error: ", errmsg(rn)))
  } else {
    ok <- same(ro$value, rn$value, tol)
    record(label, "SAME", if (ok) "PASS" else "FAIL",
      if (ok) "" else paste(utils::head(all.equal(ro$value, rn$value, tolerance = tol), 2), collapse = " / "))
  }
}

# Same, over a named list of argument lists: one summary row for all settings.
cmp_grid <- function(label, fn, arg_list, tol = 1e-8) {
  bad <- character()
  for (i in seq_along(arg_list)) {
    a <- arg_list[[i]]
    ro <- run(do.call(OLD[[fn]], a)); rn <- run(do.call(NEW[[fn]], a))
    ok <- if (ro$error && rn$error) {
      TRUE
    } else if (ro$error || rn$error) {
      FALSE
    } else {
      same(ro$value, rn$value, tol)
    }
    if (!ok) bad <- c(bad, if (is.null(names(arg_list))) as.character(i) else names(arg_list)[i])
  }
  record(sprintf("%s [%d settings]", label, length(arg_list)), "SAME",
    if (length(bad) == 0) "PASS" else "FAIL",
    if (length(bad) == 0) "" else paste("differs for:", paste(bad, collapse = "; ")))
}

# Expect a documented difference.
expect_diff <- function(label, cond, detail = "") {
  ok <- tryCatch(isTRUE(cond), error = function(e) FALSE)
  record(label, "INTENDED_DIFF", if (ok) "PASS" else "FAIL", detail)
}
skip <- function(label, why) record(label, "SKIP", "SKIP", why)


# =============================================================================
# hat_matrix.R : hat_Q (new), hat_matrix
# =============================================================================
covers("hat_Q", "hat_matrix")
set.seed(1)
T_ <- 100
D1 <- cbind(1, seq_len(T_), rnorm(T_), rnorm(T_))
D2 <- cbind(D1, D1[, 3] + D1[, 4])                          # collinear column
spk <- NEW$flags2spikes(c(5, 6, 40), T_)
D3 <- cbind(D1, NEW$dct_bases(T_, 6), spk, spk[, 1])        # DCT + spikes + duplicate
for (nm in c("D1", "D2", "D3")) {
  D <- get(nm)
  cmp(paste("hat_matrix", nm), OLD$hat_matrix(D), NEW$hat_matrix(D), tol = 1e-10)
  Qd <- NEW$hat_Q(D)
  expect_diff(paste("hat_Q orthonormal", nm),
    same(crossprod(Qd), diag(ncol(Qd)), 1e-10),
    "new helper: Q'Q = I")
}
cmp("hat_matrix NA in design (both error)", OLD$hat_matrix(cbind(1, c(NA, seq_len(T_ - 1)))), NEW$hat_matrix(cbind(1, c(NA, seq_len(T_ - 1)))))
r_old <- run(OLD$hat_matrix(matrix(0, T_, 1))); r_new <- run(NEW$hat_matrix(matrix(0, T_, 1)))
expect_diff("hat_matrix: all-zero design (rank 0) gives the zero projector in new",
  !r_new$error && same(r_new$value, matrix(0, T_, T_)),
  paste("old:", if (r_old$error) errmsg(r_old) else "returned a value"))


# =============================================================================
# nuisance_regression.R
# =============================================================================
covers("nuisance_regression")
set.seed(2)
Y_TV <- matrix(rnorm(T_ * 30), T_, 30)
for (nm in c("D1", "D2", "D3")) {
  D <- get(nm)
  cmp(paste("nuisance_regression T x V,", nm), OLD$nuisance_regression(Y_TV, D), NEW$nuisance_regression(Y_TV, D), 1e-9)
  cmp(paste("nuisance_regression V x T,", nm), OLD$nuisance_regression(t(Y_TV), D), NEW$nuisance_regression(t(Y_TV), D), 1e-9)
}
Y_sq <- matrix(rnorm(T_ * T_), T_, T_)
Y_na <- Y_TV; Y_na[3, 2] <- NA
cmp("nuisance_regression square Y (treated as T x V)", OLD$nuisance_regression(Y_sq, D1), NEW$nuisance_regression(Y_sq, D1), 1e-9)
cmp("nuisance_regression NA in Y", OLD$nuisance_regression(Y_na, D3), NEW$nuisance_regression(Y_na, D3), 1e-9)
cmp("nuisance_regression single-vector Y", OLD$nuisance_regression(Y_TV[, 1], D3), NEW$nuisance_regression(Y_TV[, 1], D3), 1e-9)
cmp("nuisance_regression incompatible dims (both error)", OLD$nuisance_regression(matrix(1, 7, 8), D1), NEW$nuisance_regression(matrix(1, 7, 8), D1))
# Residuals are orthogonal to the design (sanity check on the new version)
expect_diff("nuisance_regression residuals orthogonal to design",
  max(abs(crossprod(D1, NEW$nuisance_regression(Y_TV, D1)))) < 1e-8, "sanity")


# =============================================================================
# DCT_FFT_detrending.R : dct_bases, dct2Hz, Hz2dct, fft_detrend, temporal_filter
# =============================================================================
covers("dct_bases", "dct2Hz", "Hz2dct", "fft_detrend", "temporal_filter", "detrend")
cmp_grid("dct_bases", "dct_bases", list(
  "T=50, n=1" = list(T_ = 50, n = 1), "T=50, n=10" = list(T_ = 50, n = 10),
  "T=50, n=50" = list(T_ = 50, n = 50), "T=7, n=7" = list(T_ = 7, n = 7),
  "T=7, n=3" = list(T_ = 7, n = 3), "T=1200, n=17" = list(T_ = 1200, n = 17),
  "n > T (both error)" = list(T_ = 10, n = 11), "n = 0 (both error)" = list(T_ = 10, n = 0)
), 1e-10)
cmp("dct_bases idx= matches columns of old full matrix",
    OLD$dct_bases(50, 50)[, c(3, 10, 50)], NEW$dct_bases(50, idx = c(3, 10, 50)), 1e-10)
cmp_grid("dct2Hz", "dct2Hz", list(list(T_ = 1200, TR = .72, n = 17), list(T_ = 200, TR = 2, n = 3)))
cmp_grid("Hz2dct", "Hz2dct", list(list(T_ = 1200, TR = .72, f = .01), list(T_ = 200, TR = 2, f = .005)))

# fft_detrend: intended difference (conjugate-symmetric zeroing) ---------------
# Old zeroed bins 1..N and T-N+1..T, which is not conjugate-symmetric. After
#  Re(), the only bins that change relative to new are the +N / -N cycle bins,
#  which old halves; everything else is identical.
set.seed(3)
T2 <- 200
x <- apply(matrix(rnorm(T2 * 3), T2, 3), 2, cumsum)
Fx <- mvfft(x)
for (N in c(1, 2, 4, 10)) {
  Fo <- mvfft(OLD$fft_detrend(x, N)); Fn <- mvfft(NEW$fft_detrend(x, N))
  rm_bins <- c(seq_len(N), if (N > 1) seq(T2 - N + 2, T2))
  half <- c(N + 1, T2 - N + 1)
  other <- setdiff(seq_len(T2), c(rm_bins, half))
  expect_diff(sprintf("fft_detrend N=%d: new zeroes bins 0..N-1 (+conjugates), keeps all others as in input", N),
    all(Mod(Fn[rm_bins, ]) < 1e-8) && same(Fn[-rm_bins, ], Fx[-rm_bins, ], 1e-8), "")
  expect_diff(sprintf("fft_detrend N=%d: old differs only at the +/-N bins, which it halves", N),
    all(Mod(Fo[rm_bins, ]) < 1e-8) && same(Fo[other, ], Fx[other, ], 1e-8) && same(Fo[half, ], Fx[half, ] / 2, 1e-8), "")
}

# temporal_filter ---------------------------------------------------------------
set.seed(4)
T3 <- 1200; TR <- .72
X3 <- matrix(rnorm(T3 * 40), T3, 40) + 0.01 * seq_len(T3)
cmp("temporal_filter DCT hpf only",        OLD$temporal_filter(X3, TR, hpf = .01),            NEW$temporal_filter(X3, TR, hpf = .01), 1e-8)
cmp("temporal_filter DCT default hpf",     OLD$temporal_filter(X3, TR),                       NEW$temporal_filter(X3, TR), 1e-8)
cmp("temporal_filter DCT hpf + lpf",       OLD$temporal_filter(X3, TR, hpf = .01, lpf = .2),  NEW$temporal_filter(X3, TR, hpf = .01, lpf = .2), 1e-8)
cmp("temporal_filter DCT lpf only",        OLD$temporal_filter(X3, TR, hpf = NULL, lpf = .2), NEW$temporal_filter(X3, TR, hpf = NULL, lpf = .2), 1e-8)
cmp("temporal_filter returns bases (X = T)", OLD$temporal_filter(T3, TR, hpf = .01),          NEW$temporal_filter(T3, TR, hpf = .01), 1e-10)
cmp("temporal_filter bases hpf + lpf (X = T)", OLD$temporal_filter(T3, TR, hpf = .01, lpf = .2), NEW$temporal_filter(T3, TR, hpf = .01, lpf = .2), 1e-10)
cmp("temporal_filter hpf = lpf = NULL (matrix)", OLD$temporal_filter(X3, TR, hpf = NULL, lpf = NULL), NEW$temporal_filter(X3, TR, hpf = NULL, lpf = NULL))
cmp("detrend wrapper", OLD$detrend(X3, TR), NEW$detrend(X3, TR), 1e-8)
# Grid over scan length, TR, and hpf/lpf (data, and bases when X = T).
fmt <- function(z) if (is.null(z)) "NULL" else format(z)
Xg <- list("150" = matrix(rnorm(150 * 5), 150, 5), "600" = matrix(rnorm(600 * 5), 600, 5))
tf_data <- list(); tf_bases <- list()
for (Tg in c(150, 600)) for (TRg in c(.72, 2)) for (hp in list(NULL, .008, .02)) for (lp in list(NULL, .1)) {
  if (is.null(hp) && is.null(lp)) next
  nm_ <- sprintf("T=%d, TR=%g, hpf=%s, lpf=%s", Tg, TRg, fmt(hp), fmt(lp))
  tf_data[[nm_]]  <- list(X = Xg[[as.character(Tg)]], TR = TRg, hpf = hp, lpf = lp)
  tf_bases[[nm_]] <- list(X = Tg, TR = TRg, hpf = hp, lpf = lp)
}
cmp_grid("temporal_filter DCT (data)",  "temporal_filter", tf_data,  1e-8)
cmp_grid("temporal_filter DCT (bases)", "temporal_filter", tf_bases, 1e-10)

# Intended differences
Xs <- matrix(rnorm(20 * 5), 20, 5)   # short scan
short <- lapply(c(10, 20, 30), function(n_) matrix(rnorm(n_ * 3), n_, 3))
r_sh <- lapply(short, function(Xq) list(o = run(OLD$temporal_filter(Xq, TR)), n = run(NEW$temporal_filter(Xq, TR))))
expect_diff("temporal_filter: hpf rounds to 0 bases at T = 10, 20, 30 (old errors, new = centered X)",
  all(vapply(seq_along(short), function(i)
    r_sh[[i]]$o$error && !r_sh[[i]]$n$error && same(r_sh[[i]]$n$value, NEW$colCenter(short[[i]]), 1e-10), logical(1))),
  paste("old error:", errmsg(r_sh[[1]]$o)))
cmp("temporal_filter: 0 hpf bases + lpf (unchanged)", OLD$temporal_filter(Xs, TR, hpf = .008, lpf = .15), NEW$temporal_filter(Xs, TR, hpf = .008, lpf = .15), 1e-8)
r_old <- run(OLD$temporal_filter(Xs, TR, hpf = .008, lpf = .009))
r_new <- run(NEW$temporal_filter(Xs, TR, hpf = .008, lpf = .009))
expect_diff("temporal_filter: lpf rounds to 0 bases with hpf also 0 (new stops 'too low')",
  !r_old$error && r_new$error, paste("new error:", errmsg(r_new)))
r_old <- run(OLD$temporal_filter(5, TR, hpf = NULL, lpf = NULL))
r_new <- run(NEW$temporal_filter(5, TR, hpf = NULL, lpf = NULL))
expect_diff("temporal_filter: integer X with hpf = lpf = NULL (old errors, new = T x 0)",
  r_old$error && !r_new$error && identical(dim(r_new$value), c(5L, 0L)), "")
# FFT path: differs from old only through fft_detrend.
N_fft <- round(.01 * TR * T3)
expect_diff("temporal_filter FFT: = fft_detrend(N) in both, so it differs from old only via fft_detrend",
  same(NEW$temporal_filter(X3, TR, hpf = .01, method = "FFT"), NEW$fft_detrend(X3, N_fft), 1e-10) &&
    same(OLD$temporal_filter(X3, TR, hpf = .01, method = "FFT"), OLD$fft_detrend(X3, N_fft), 1e-10) &&
    !same(OLD$temporal_filter(X3, TR, hpf = .01, method = "FFT"), NEW$temporal_filter(X3, TR, hpf = .01, method = "FFT"), 1e-6), "")
cmp("temporal_filter FFT, N < 1 (just centers)", OLD$temporal_filter(X3, TR, hpf = .0001, method = "FFT"), NEW$temporal_filter(X3, TR, hpf = .0001, method = "FFT"))
cmp("temporal_filter FFT with lpf (both error)", OLD$temporal_filter(X3, TR, hpf = .01, lpf = .2, method = "FFT"), NEW$temporal_filter(X3, TR, hpf = .01, lpf = .2, method = "FFT"))
cmp("temporal_filter FFT with integer X (both error)", OLD$temporal_filter(T3, TR, hpf = .01, method = "FFT"), NEW$temporal_filter(T3, TR, hpf = .01, method = "FFT"))


# =============================================================================
# flags2spikes.R
# =============================================================================
covers("flags2spikes")
lg <- rep(FALSE, 100); lg[c(2, 7, 99)] <- TRUE
cmp_grid("flags2spikes (sorted unique flags)", "flags2spikes", list(
  "3 flags" = list(flags = c(3, 10, 50), n_time = 100),
  "both ends" = list(flags = c(1, 100), n_time = 100),
  "every volume" = list(flags = 1:20, n_time = 20),
  "single flag" = list(flags = 7, n_time = 7),
  "logical" = list(flags = lg, n_time = 100),
  "all-TRUE logical" = list(flags = rep(TRUE, 10), n_time = 10),
  "out of range (both error)" = list(flags = c(1, 101), n_time = 100),
  "non-integer (both error)" = list(flags = 1.5, n_time = 100),
  "NA flag (both error)" = list(flags = c(1, NA), n_time = 100)
))
# Unsorted / duplicated: same projection, but old keeps the input column order
#  (and duplicates) while new sorts and de-duplicates.
cmp("flags2spikes unsorted flags: same projection",
    OLD$hat_matrix(OLD$flags2spikes(c(50, 3, 10), 100)),
    NEW$hat_matrix(NEW$flags2spikes(c(50, 3, 10), 100)), 1e-10)
expect_diff("flags2spikes: unsorted flags -> old keeps input column order, new sorts (same columns)",
  !same(OLD$flags2spikes(c(50, 3, 10), 100), NEW$flags2spikes(c(50, 3, 10), 100)) &&
    same(OLD$flags2spikes(c(50, 3, 10), 100)[, order(c(50, 3, 10))], NEW$flags2spikes(c(50, 3, 10), 100)), "")
expect_diff("flags2spikes: duplicate flags dedup'd in new (same projection)",
  ncol(OLD$flags2spikes(c(5, 5, 9), 100)) == 3 && ncol(NEW$flags2spikes(c(5, 5, 9), 100)) == 2 &&
    same(OLD$hat_matrix(OLD$flags2spikes(c(5, 5, 9), 100)), NEW$hat_matrix(NEW$flags2spikes(c(5, 5, 9), 100)), 1e-10), "")
r_old <- run(OLD$flags2spikes(integer(0), 100)); r_new <- run(NEW$flags2spikes(integer(0), 100))
expect_diff("flags2spikes: empty flags (old warns, new silent T x 0)",
  length(r_old$warnings) > 0 && length(r_new$warnings) == 0 && identical(dim(r_new$value), c(100L, 0L)),
  paste("old warnings:", length(r_old$warnings)))
expect_diff("flags2spikes: logical of wrong length (old accepts, new errors)",
  !run(OLD$flags2spikes(c(TRUE, FALSE, TRUE), 100))$error && run(NEW$flags2spikes(c(TRUE, FALSE, TRUE), 100))$error, "")


# =============================================================================
# PCA.R  (Q is supplied, so pesel is not required)
# =============================================================================
covers("PCA")
set.seed(5)
V <- 400; Tn <- 30; Q <- 5
Xp <- matrix(rnorm(V * Tn), V, Tn) %*% diag(seq(3, 1, length.out = Tn))  # distinct eigenvalues
cn <- function(M) sqrt(colSums(M^2))
o <- run(OLD$PCA(Xp, Q = Q)); n <- run(NEW$PCA(Xp, Q = Q))
if (!o$error && !n$error) {
  record("PCA d (eigenvalues)", "SAME", if (same(o$value$d, n$value$d, 1e-8)) "PASS" else "FAIL")
  record("PCA u (up to sign)", "SAME",
    if (same(abs(diag(crossprod(o$value$u, n$value$u))), rep(1, Q), 1e-6)) "PASS" else "FAIL")
  record("PCA v is NULL when Vdim = 0", "SAME",
    if (is.null(o$value$v) && is.null(n$value$v)) "PASS" else "FAIL")
} else {
  record("PCA basic", "SAME", "FAIL", paste("old:", errmsg(o), "| new:", errmsg(n)))
}
Qs <- c(1, 2, 10, 20)
pca_ok <- vapply(Qs, function(q) {
  oq <- run(OLD$PCA(Xp, Q = q)); nq <- run(NEW$PCA(Xp, Q = q))
  !oq$error && !nq$error && same(oq$value$d, nq$value$d, 1e-8) &&
    same(abs(diag(crossprod(oq$value$u, nq$value$u))), rep(1, q), 1e-6)
}, logical(1))
record(sprintf("PCA d and u (up to sign) [%d settings: Q = %s]", length(Qs), paste(Qs, collapse = ", ")),
       "SAME", if (all(pca_ok)) "PASS" else "FAIL", if (all(pca_ok)) "" else paste("differs for Q =", paste(Qs[!pca_ok], collapse = ", ")))
Xc <- scale(Xp, scale = FALSE)
cmp("PCA center = FALSE on pre-centered data (d)",
    OLD$PCA(Xc, center = FALSE, Q = Q)$d, NEW$PCA(Xc, center = FALSE, Q = Q)$d, 1e-8)
cmp("PCA center = FALSE on uncentered data (both error)", OLD$PCA(Xp + 5, center = FALSE, Q = Q), NEW$PCA(Xp + 5, center = FALSE, Q = Q))
cmp("PCA data.frame input (both error)", OLD$PCA(as.data.frame(Xp), Q = Q), NEW$PCA(as.data.frame(Xp), Q = Q))
Xw <- matrix(rnorm(10 * Tn), 10, Tn)
cmp("PCA more time points than voxels: d", OLD$PCA(Xw, Q = 3)$d, NEW$PCA(Xw, Q = 3)$d, 1e-8)
cmp("PCA Vdim > Q: capped to Q (dims of v)", dim(OLD$PCA(Xp, Q = 3, Vdim = 5)$v), dim(NEW$PCA(Xp, Q = 3, Vdim = 5)$v))

# Directions: same direction up to sign, but old is not unit norm.
vd <- lapply(list(2, 3, 5, "Q"), function(v_) list(
  o = run(OLD$PCA(Xp, Q = Q, Vdim = v_)), n = run(NEW$PCA(Xp, Q = Q, Vdim = v_)),
  k = if (identical(v_, "Q")) Q else v_))
vd_all <- function(f) all(vapply(vd, function(z) !z$o$error && !z$n$error && f(z), logical(1)))
expect_diff("PCA Vdim = 2, 3, 5, 'Q': directions parallel (up to sign)",
  vd_all(function(z) same(abs(colSums(z$o$value$v * z$n$value$v)) / (cn(z$o$value$v) * cn(z$n$value$v)), rep(1, z$k), 1e-6)), "")
expect_diff("PCA Vdim = 2, 3, 5, 'Q': new directions are unit norm",
  vd_all(function(z) same(cn(z$n$value$v), rep(1, z$k), 1e-8)), "")
expect_diff("PCA Vdim = 2, 3, 5, 'Q': old norms were sqrt((nV-1)/d), not 1",
  vd_all(function(z) same(cn(z$o$value$v), sqrt((V - 1) / z$n$value$d[seq_len(z$k)]), 1e-6)), "")
o <- run(OLD$PCA(Xp, Q = Q, Vdim = 1)); n <- run(NEW$PCA(Xp, Q = Q, Vdim = 1))
expect_diff("PCA Vdim = 1: old errors / wrong shape (diag(scalar) bug), new is V x 1 unit norm",
  (o$error || !identical(dim(o$value$v), c(as.integer(V), 1L))) && !n$error &&
    identical(dim(n$value$v), c(as.integer(V), 1L)) && same(cn(n$value$v), 1, 1e-8),
  paste("old:", errmsg(o)))
r_old <- run(OLD$PCA(Xp, Q = Tn + 3)); r_new <- run(NEW$PCA(Xp, Q = Tn + 3))
expect_diff("PCA Q > number of columns: new caps Q at ncol(X) with a warning",
  !r_new$error && ncol(r_new$value$u) == Tn && length(r_new$warnings) > 0,
  paste("old:", if (r_old$error) errmsg(r_old) else paste0("returned u with ", ncol(r_old$value$u), " columns")))
if (has("pesel")) {
  cmp("PCA with pesel-selected Q (d)", OLD$PCA(Xp)$d, NEW$PCA(Xp)$d, 1e-8)
} else skip("PCA with pesel-selected Q", "pesel not installed")


# =============================================================================
# despike_3D.R
# =============================================================================
# A series whose trend lies in both the old and the new reference functions;
#  bounded noise keeps every |s| below c1, so nothing is edited unless spiked.
mk_trend <- function(n_) { tt <- seq_len(n_); 3 * cos(2 * pi * tt / n_) + sin(2 * pi * tt / n_) }
mk_y <- function(n_) runif(n_, -1, 1) + mk_trend(n_)

# New-only checks (no fda / quantreg needed) ------------------------------------
covers("despike_3D", "despike_3D.fit_NEW")
set.seed(6)
y6 <- mk_y(600); sp6 <- c(100, 300, 500)
y6s <- y6; y6s[sp6] <- y6s[sp6] + c(50, -50, 50)
o6 <- NEW$despike_3D(y6s)
expect_diff("despike_3D method = 'auto': least-squares NEW fit for > 500 points (= method = 'NEW')",
  same(o6, NEW$despike_3D(y6s, method = "NEW"), 1e-10), "")
expect_diff("despike_3D NEW: compresses spikes of both signs, leaves every other point untouched",
  same(o6[-sp6], y6s[-sp6], 1e-10) && all(abs(o6[sp6] - mk_trend(600)[sp6]) < 3.5), "")
y3 <- mk_y(300); y3[c(40, 90)] <- y3[c(40, 90)] + c(50, -50)
expect_diff("despike_3D localedit: spike = average of its non-spike neighbours; others untouched",
  same(NEW$despike_3D(y3, method = "NEW", localedit = TRUE)[c(40, 90)], c(mean(y3[c(39, 41)]), mean(y3[c(89, 91)])), 1e-10), "")
expect_diff("despike_3D: < 15 time points, c1 < 1, or c2 < c1 + 0.5 now error",
  run(NEW$despike_3D(rnorm(10)))$error && run(NEW$despike_3D(rnorm(50), c1 = .5))$error &&
    run(NEW$despike_3D(rnorm(50), c1 = 3, c2 = 3.2))$error, "")

# Old vs new -----------------------------------------------------------------------
if (has("fda")) {
  # Reference functions: old Fourier basis is contained in the new one.
  span_ok <- vapply(c(40, 100, 300, 600), function(nT) {
    old_basis <- fda::getbasismatrix(seq(nT), fda::create.fourier.basis(
      rangeval = c(0, nT), nbasis = 2 * round((nT / 30 / 2)) + 1))
    new_basis <- NEW$despike_3D.basis(nT)
    ncol(new_basis) == 3 + 2 * round(nT / 30) && qr(cbind(new_basis, old_basis))$rank == qr(new_basis)$rank
  }, logical(1))
  covers("despike_3D.basis")
  expect_diff("despike basis (nT = 40, 100, 300, 600): new = const + linear + quadratic + round(nT/30) sin/cos pairs, containing the old Fourier basis",
    all(span_ok), "")
} else skip("despike basis vs fda", "fda not installed")

if (has("fda") && has("quantreg")) {
  covers("despike_3D", "despike_3D.qreg", "despike_3D.interpolate")
  set.seed(6)

  # (1) No spikes: old/new differ only at marginal points.
  for (n_ in c(60, 150, 300)) {
    yy <- mk_y(n_)
    o <- as.numeric(OLD$despike_3D(yy)); n <- as.numeric(NEW$despike_3D(yy))
    d <- abs(o - n)
    expect_diff(sprintf("despike_3D no spike, n = %d: old/new differ only at marginal points", n_),
                mean(d > 1e-8) < 0.05 && max(d) < 0.1,
                sprintf("%d points differ, max diff %.3g", sum(d > 1e-8), max(d)))
  }

  # (2) Positive spikes
  y <- mk_y(300); sp <- c(40, 90, 150, 200, 260)
  yp <- y; yp[sp] <- yp[sp] + 50
  o <- as.numeric(OLD$despike_3D(yp)); n <- as.numeric(NEW$despike_3D(yp))

  expect_diff("despike_3D positive spikes: old and new each edit exactly the spike points",
    same(o[-sp], yp[-sp], 1e-2) && same(n[-sp], yp[-sp], 1e-2) && all(o[sp] < yp[sp] - 40) && all(n[sp] < yp[sp] - 40), "")

  # (3) Sigma: a saturated edit sits ~c2*sigma above the fit, so the ratio of edit heights is
  #     sigma_new / sigma_old = sqrt(pi/2) * median|r| / (1.4826 * MAD(r)) ~ 0.85.
  ratio <- mean((n[sp] - mk_trend(300)[sp]) / (o[sp] - mk_trend(300)[sp]))
  expect_diff("despike_3D sigma: edited spikes sit lower with new sigma = sqrt(pi/2)*median|r| (~0.85x the old 1.4826*MAD)",
    ratio > .75 && ratio < .95, sprintf("mean ratio of edit heights new/old = %.3f (expect ~%.3f)", ratio, sqrt(pi / 2) / 1.4826))

  # (4) Negative spike: old leaves it, new compresses it; no other point differs.
  yn <- y; yn[80] <- yn[80] - 50
  o <- as.numeric(OLD$despike_3D(yn)); n <- as.numeric(NEW$despike_3D(yn))
  expect_diff("despike_3D negative spike: old leaves it, new compresses it (every other point untouched in both)",
    same(o[80], yn[80], 1e-8) && n[80] > yn[80] + 40 &&
      same(o[-80], yn[-80], 1e-2) && same(n[-80], yn[-80], 1e-2), "")

  # (5) Default method: L1 (OLD) up to 500 time points.
  expect_diff("despike_3D method = 'auto': L1 fit (= method = 'OLD') for <= 500 points",
    same(NEW$despike_3D(yp), NEW$despike_3D(yp, method = "OLD"), 1e-10), "")

  # (6) Output and degenerate input.
  oo <- run(OLD$despike_3D(yp)); nn <- run(NEW$despike_3D(yp))
  expect_diff("despike_3D output: old kept the names of the fitted values, new is a plain numeric vector",
    !is.null(names(oo$value)) && is.null(names(nn$value)), "")
  yc <- rep(1, 50)
  oo <- run(OLD$despike_3D(yc)); nn <- run(NEW$despike_3D(yc))
  expect_diff("despike_3D constant series (MAD = 0): old NaN/NA or error, new returns input",
    (oo$error || anyNA(oo$value)) && !nn$error && same(as.numeric(nn$value), yc, 1e-10), paste("old:", errmsg(oo)))
} else skip("despike_3D old-vs-new comparisons", "fda and/or quantreg not installed")


# =============================================================================
# utils_mat.R : rowMedians2, unvec_mat, pct_sig
# =============================================================================
covers("rowMedians2", "unvec_mat", "pct_sig")
set.seed(7)
xm <- matrix(rnorm(50 * 12), 50, 12); xm[3, 4] <- NA
robust <- has("robustbase")
r_old <- run(OLD$rowMedians2(xm, na.rm = TRUE)); r_new <- run(NEW$rowMedians2(xm, na.rm = TRUE))
ref <- apply(xm, 1, median, na.rm = TRUE)
expect_diff("rowMedians2: new matches apply(median) with and without na.rm",
  !r_new$error && same(as.numeric(r_new$value), ref) &&
    same(as.numeric(NEW$rowMedians2(xm)), apply(xm, 1, median)), "")
expect_diff("rowMedians2: old inverted condition (errors w/o robustbase, else apply)",
  if (robust) (!r_old$error && same(as.numeric(r_old$value), ref)) else r_old$error,
  paste("robustbase installed:", robust, "| old:", errmsg(r_old)))

xv <- cbind(seq(3), seq(2, 4), seq(3, 5)); mk <- matrix(c(rep(TRUE, 3), FALSE), ncol = 2)
mk2 <- matrix(rbinom(20 * 15, 1, .6) == 1, 20, 15, dimnames = list(letters[1:20], LETTERS[1:15]))
xr <- matrix(rnorm(sum(mk2) * 4), sum(mk2), 4)
cmp_grid("unvec_mat", "unvec_mat", list(
  "docs example (fill 0)" = list(x = xv, mask = mk, fill_value = 0),
  "fill NA" = list(x = xv, mask = mk),
  "numeric 0/1 mask" = list(x = xv, mask = 1 * mk, fill_value = 0),
  "random mask, dimnames, fill -1" = list(x = xr, mask = mk2, fill_value = -1),
  "random mask, fill NA" = list(x = xr, mask = mk2),
  "integer data, fill NA" = list(x = matrix(seq_len(sum(mk2)), ncol = 1), mask = mk2),
  "all-TRUE mask, one column" = list(x = matrix(rnorm(6), 6, 1), mask = matrix(TRUE, 2, 3)),
  "wrong number of rows (both error)" = list(x = xv[1:2, ], mask = mk)
))
vec <- seq_len(sum(mk2)) * 1.5
r_old <- run(OLD$unvec_mat(vec, mk2)); r_new <- run(NEW$unvec_mat(vec, mk2))
expect_diff("unvec_mat plain vector: old errors, new = old on as.matrix(vec)",
  r_old$error && !r_new$error && same(r_new$value, OLD$unvec_mat(as.matrix(vec), mk2)), paste("old:", errmsg(r_old)))

xp <- matrix(runif(40 * 6, 50, 150), 40, 6, dimnames = list(NULL, paste0("v", 1:6)))
cmp_grid("pct_sig", "pct_sig", list(
  "median/column" = list(X = xp), "mean/column" = list(X = xp, center = mean),
  "median/all" = list(X = xp, by = "all"), "mean/all" = list(X = xp, center = mean, by = "all"),
  "max/column" = list(X = xp, center = max), "no dimnames" = list(X = unname(xp))
))


# =============================================================================
# scale_matrix.R : scale_med, scale_timeseries, scale_design_mat
# (OLD$rowMedians2 is replaced by the correct apply version so scale_med is
#  compared in isolation from the rowMedians2 bug.)
# =============================================================================
covers("scale_med", "scale_timeseries", "scale_design_mat")
OLD$rowMedians2 <- function(x, na.rm = FALSE, ...) apply(x, 1, median, na.rm = na.rm, ...)

set.seed(8)
M <- matrix(rnorm(60 * 8, 5, 3), 60, 8)
Mna <- M; Mna[4, 2] <- NA
Mc <- M; Mc[, 3] <- 2; Mc[, 6] <- 7
cmp_grid("scale_med", "scale_med", list(
  "columns" = list(mat = M), "with NA entries" = list(mat = Mna),
  "2 const cols, drop" = list(mat = Mc), "2 const cols, keep as NA" = list(mat = Mc, drop_const = FALSE),
  "TOL = 1e-2" = list(mat = M, TOL = 1e-2), "all const (both error)" = list(mat = matrix(1, 10, 3))
))
# Intended differences
r_old <- run(OLD$scale_med(M, doRows = TRUE)); r_new <- run(NEW$scale_med(M, doRows = TRUE))
expect_diff("scale_med doRows: old returns the transpose (orientation bug); new keeps input orientation, same values",
  !r_old$error && !r_new$error && identical(dim(r_new$value), dim(M)) && same(t(r_old$value), r_new$value), "")
M1 <- cbind(2, 2, rnorm(50))          # exactly one non-constant column
r_old <- run(OLD$scale_med(M1)); r_new <- run(NEW$scale_med(M1))
expect_diff("scale_med: single non-constant col -> old is 1 x T (wrong), new is T x 1; values equal",
  !r_old$error && !r_new$error && identical(dim(r_old$value), c(1L, 50L)) &&
    identical(dim(r_new$value), c(50L, 1L)) && same(as.numeric(r_old$value), as.numeric(r_new$value)), "")
Mn <- M; Mn[, 5] <- NA
r_old <- run(OLD$scale_med(Mn)); r_new <- run(NEW$scale_med(Mn))
expect_diff("scale_med: all-NA column -> old errors, new treats as constant and drops it",
  r_old$error && !r_new$error && ncol(r_new$value) == 7, paste("old:", errmsg(r_old)))
expect_diff("scale_med warning text: with drop_const = FALSE new says 'set to NA' (old said 'removed')",
  any(grepl("removed", run(OLD$scale_med(Mc, drop_const = FALSE))$warnings)) &&
    any(grepl("set to NA", run(NEW$scale_med(Mc, drop_const = FALSE))$warnings)), "")

# scale_timeseries (V > T so the default transpose logic does not fire, except
#  for the T x V and "V < T" cases noted below)
set.seed(8)
Bm <- matrix(rnorm(100 * 60, 100, 5), 100, 60)     # mean scaling under "auto"
Bs <- matrix(rnorm(100 * 60), 100, 60)             # sd scaling under "auto"
Bna <- Bs; Bna[3, 7] <- NA; Bna[10, 20:25] <- NA
Bmna <- Bm; Bmna[2, 3] <- NA
Bp <- matrix(rnorm(40 * 120, 50, 2), 40, 120)      # parcellated-style: V < T
Bz <- Bs; Bz[5, ] <- 0
Bsn <- Bs; Bsn[4, ] <- NA
Blow <- matrix(rnorm(100 * 60, .05, .001), 100, 60)
cmp_grid("scale_timeseries", "scale_timeseries", list(
  "auto -> mean" = list(BOLD = Bm), "auto -> sd" = list(BOLD = Bs),
  "mean" = list(BOLD = Bm, scale = "mean"), "mean with NAs" = list(BOLD = Bmna, scale = "mean"),
  "sd" = list(BOLD = Bm, scale = "sd"), "sd with NAs" = list(BOLD = Bna, scale = "sd"),
  "none" = list(BOLD = Bm, scale = "none"),
  "T x V input (auto-transposed)" = list(BOLD = t(Bm)),
  "V < T, transpose = FALSE" = list(BOLD = Bp, transpose = FALSE),
  "zero-variance row, sd (both error)" = list(BOLD = Bz, scale = "sd"),
  "all-NA row, sd (both error)" = list(BOLD = Bsn, scale = "sd"),
  "means < 0.1, mean (both error)" = list(BOLD = Blow, scale = "mean")
), 1e-10)

# scale_design_mat
set.seed(9)
Xd <- cbind(1, abs(rnorm(30)) + 1, runif(30, 1, 2), seq_len(30))     # positive max-magnitude
cmp_grid("scale_design_mat (positive-dominated)", "scale_design_mat", list(
  "columns" = list(x = Xd), "rows" = list(x = Xd, doRows = TRUE),
  "t(X), rows" = list(x = t(Xd), doRows = TRUE)))
Xm <- matrix(rnorm(30 * 12), 30, 12)                                 # mixed-sign columns
ratio <- apply(Xm, 2, max) / apply(abs(Xm), 2, max)
expect_diff("scale_design_mat scope: new = old * max/max|.| per column; identical where the largest-magnitude value is positive",
  any(ratio == 1) && any(ratio < 1) &&
    same(NEW$scale_design_mat(Xm), OLD$scale_design_mat(Xm) * rep(ratio, each = 30), 1e-10) &&
    same(NEW$scale_design_mat(Xm)[, ratio == 1], OLD$scale_design_mat(Xm)[, ratio == 1], 1e-10), "")
Xn <- cbind(1, -abs(rnorm(30)) - 5)                                  # negative-dominated column
expect_diff("scale_design_mat: negative-dominated column has opposite sign in old vs new",
  cor(OLD$scale_design_mat(Xn)[, 2], NEW$scale_design_mat(Xn)[, 2]) < -0.99, "")
Xz <- cbind(Xd, 0)
r_old <- run(OLD$scale_design_mat(Xz)); r_new <- run(NEW$scale_design_mat(Xz))
expect_diff("scale_design_mat: all-zero column -> old NaN, new zeros",
  !r_old$error && any(is.nan(r_old$value)) && !r_new$error && all(r_new$value[, 5] == 0), "")


# =============================================================================
# mask_BOLD.R / utils_misc.R : mask_BOLD, make_mask
# =============================================================================
if (has("matrixStats")) {
  covers("mask_BOLD", "make_mask")
  set.seed(10)
  B <- matrix(rnorm(50 * 40), 50, 40)
  B[3, 5] <- NA; B[7, ] <- 1; B[9, ] <- rnorm(40, -100)
  Bnan <- B; Bnan[12, 3] <- NaN
  Bi <- B; Bi[11, 2] <- Inf
  Bi2 <- B; Bi2[11, 2] <- -Inf
  Ball <- B; Ball[, 1] <- NA
  cmp_grid("mask_BOLD", "mask_BOLD", list(
    "defaults" = list(BOLD = B, verbose = FALSE),
    "meanTol = -50" = list(BOLD = B, meanTol = -50, verbose = FALSE),
    "varTol = 0.5" = list(BOLD = B, varTol = 0.5, verbose = FALSE),
    "varTol = 0" = list(BOLD = B, varTol = 0, verbose = FALSE),
    "meanTol = -50, varTol = 0.5" = list(BOLD = B, meanTol = -50, varTol = 0.5, verbose = FALSE),
    "no NA rows" = list(BOLD = B[-3, ], verbose = FALSE),
    "every row has NA" = list(BOLD = Ball, verbose = FALSE),
    "NaN entry" = list(BOLD = Bnan, verbose = FALSE),
    "+Inf row, varTol = 0" = list(BOLD = Bi, varTol = 0, verbose = FALSE),
    "-Inf row, meanTol = -50" = list(BOLD = Bi2, meanTol = -50, verbose = FALSE)
  ))
  if (is.function(NEW$make_mask)) {
    cmp("make_mask (old duplicate vs new wrapper)", OLD$make_mask(B, verbose = FALSE), NEW$make_mask(B, verbose = FALSE))
    cmp("make_mask == mask_BOLD in new", NEW$mask_BOLD(B, verbose = FALSE), NEW$make_mask(B, verbose = FALSE))
  } else {
    # make_mask was removed in new; its job is done by mask_BOLD.
    cmp("make_mask removed in new: old make_mask == new mask_BOLD", OLD$make_mask(B, verbose = FALSE), NEW$mask_BOLD(B, verbose = FALSE))
  }
  # Infinite value -> NaN variance: old leaves the row TRUE, new masks it (only that row differs).
  expect_diff("mask_BOLD: row containing +Inf (NaN variance) is masked in new, kept in old; no other row differs",
    isTRUE(OLD$mask_BOLD(Bi, verbose = FALSE)[11]) && !NEW$mask_BOLD(Bi, verbose = FALSE)[11] &&
      identical(which(OLD$mask_BOLD(Bi, verbose = FALSE) != NEW$mask_BOLD(Bi, verbose = FALSE)), 11L), "")
  # Verbose wording was normalized; the counts reported are unchanged.
  nums <- function(f, ...) {
    o <- capture.output(invisible(f(B, ...)))
    as.numeric(unlist(regmatches(o, gregexpr("[0-9]+", o))))
  }
  cmp("mask_BOLD verbose: counts reported are unchanged (defaults)", nums(OLD$mask_BOLD), nums(NEW$mask_BOLD))
  cmp("mask_BOLD verbose: counts reported are unchanged (meanTol = -50)", nums(OLD$mask_BOLD, meanTol = -50), nums(NEW$mask_BOLD, meanTol = -50))
  cat("\n-- verbose output (eyeball only; wording intentionally normalized) --\nOLD: ")
  invisible(OLD$mask_BOLD(B)); cat("\nNEW: "); invisible(NEW$mask_BOLD(B)); cat("\n\n")
} else skip("mask_BOLD / make_mask", "matrixStats not installed")


# =============================================================================
# utils_checkValue.R : match_exactly, is_1, is_posNum, is_integer
# =============================================================================
covers("match_exactly", "is_1", "is_posNum", "is_integer")
cmp_grid("match_exactly (matching input)", "match_exactly", list(
  "character" = list(user = c("a", "b"), expected = c("a", "b"), fail_action = "nothing"),
  "numeric" = list(user = c(1, 2, 3), expected = c(1, 2, 3), fail_action = "nothing"),
  "length zero" = list(user = character(0), expected = character(0), fail_action = "nothing"),
  "mismatch + stop (both error)" = list(user = "a", expected = "b", fail_action = "stop")
))
for (case in list(list(c("a", "b"), c("a", "c")), list(c("a"), c("a", "b")), list(c("b", "a"), c("a", "b")))) {
  ro <- run(OLD$match_exactly(case[[1]], case[[2]], "nothing"))
  rn <- run(NEW$match_exactly(case[[1]], case[[2]], "nothing"))
  expect_diff(sprintf("match_exactly mismatch (%s vs %s): old returns TRUE (bug), new FALSE",
                      paste(case[[1]], collapse = ","), paste(case[[2]], collapse = ",")),
    isTRUE(ro$value) && identical(rn$value, FALSE), "")
}
cmp("match_exactly mismatch: message text unchanged", msgs_of(OLD$match_exactly("a", "b")), msgs_of(NEW$match_exactly("a", "b")))
expect_diff("match_exactly mismatch + 'warning' emits a warning in both",
  length(run(OLD$match_exactly("a", "b", "warning"))$warnings) == 1 &&
    length(run(NEW$match_exactly("a", "b", "warning"))$warnings) == 1, "")

# is_1 / is_posNum / is_integer: same on every non-NA input (NA/NaN/Inf below).
ok_inputs <- list(3, 0, -1, 2.5, 1e6, "a", c(1, 2), TRUE, FALSE, NULL, 1L, numeric(0))
a_is1 <- list(); a_pos <- list(); a_int <- list()
for (i in seq_along(ok_inputs)) {
  x <- ok_inputs[[i]]
  for (dt in c("numeric", "logical", "character")) a_is1[[sprintf("#%d, %s", i, dt)]] <- list(x = x, dtype = dt)
  a_pos[[sprintf("#%d", i)]]            <- list(x = x)
  a_pos[[sprintf("#%d, zero_ok", i)]]   <- list(x = x, zero_ok = TRUE)
  a_int[[sprintf("#%d", i)]]            <- list(x = x)
  a_int[[sprintf("#%d, nneg", i)]]      <- list(x = x, nneg = TRUE)
}
cmp_grid("is_1 (non-NA inputs)", "is_1", a_is1)
cmp_grid("is_posNum (non-NA inputs)", "is_posNum", a_pos)
cmp_grid("is_integer (non-NA inputs)", "is_integer", a_int)
expect_diff("is_1(NA, 'logical'): old TRUE, new FALSE", OLD$is_1(NA, "logical") && !NEW$is_1(NA, "logical"), "")
expect_diff("is_1(NA_real_): old TRUE, new FALSE",      OLD$is_1(NA_real_) && !NEW$is_1(NA_real_), "")
expect_diff("is_1(NaN): old TRUE, new FALSE",           OLD$is_1(NaN) && !NEW$is_1(NaN), "")
expect_diff("is_posNum(NA_real_): old NA, new FALSE",   is.na(OLD$is_posNum(NA_real_)) && identical(NEW$is_posNum(NA_real_), FALSE), "")
expect_diff("is_integer(NA_real_): old errors in if(), new FALSE",
  run(OLD$is_integer(NA_real_))$error && identical(NEW$is_integer(NA_real_), FALSE), "")
expect_diff("is_integer(Inf): old errors, new FALSE",
  run(OLD$is_integer(Inf))$error && identical(NEW$is_integer(Inf), FALSE), "")


# =============================================================================
# utils_misc.R : get_FORMAT, check_req_ifti_pkg
# =============================================================================
covers("get_FORMAT", "check_req_ifti_pkg")
cmp_grid("get_FORMAT", "get_FORMAT", lapply(
  setNames(c("CIFTI", "xifti", "GIFTI", "gifti", "GIFTI2", "gifti2", "NIFTI", "nifti", "RDS", "data"),
           c("CIFTI", "xifti", "GIFTI", "gifti", "GIFTI2", "gifti2", "NIFTI", "nifti", "RDS", "data")),
  function(f) list(format = f)))
expect_diff("get_FORMAT('bogus'): old NULL, new error",
  is.null(OLD$get_FORMAT("bogus")) && run(NEW$get_FORMAT("bogus"))$error, "")
for (f in c("CIFTI", "GIFTI", "GIFTI2", "NIFTI", "MATRIX")) {
  ro <- run(OLD$check_req_ifti_pkg(f)); rn <- run(NEW$check_req_ifti_pkg(f))
  if (!ro$error && !rn$error) {
    record(paste("check_req_ifti_pkg", f), "SAME", if (same(ro$value, rn$value)) "PASS" else "FAIL")
  } else if (ro$error && rn$error) {
    msg_ok <- if (f %in% c("GIFTI", "GIFTI2")) !grepl("NIFTI|CIFTI data", rn$value$msg) || grepl("GIFTI data", rn$value$msg) else identical(ro$value$msg, rn$value$msg)
    record(paste("check_req_ifti_pkg", f, "(pkg missing: message)"), "INTENDED_DIFF", if (msg_ok) "PASS" else "FAIL",
           paste("old:", ro$value$msg, "| new:", rn$value$msg))
  } else {
    record(paste("check_req_ifti_pkg", f), "SAME", "FAIL", paste("old:", errmsg(ro), "| new:", errmsg(rn)))
  }
}


# =============================================================================
# bandstop_filter.R
# =============================================================================
if (has("gsignal")) {
  covers("bandstop_filter")
  set.seed(11)
  Xb <- matrix(rnorm(200 * 3), 200, 3)
  v <- rnorm(100)
  cmp_grid("bandstop_filter", "bandstop_filter", list(
    "matrix, TR=.72, .31-.43" = list(X = Xb, TR = .72, f1 = .31, f2 = .43),
    "TR=.8, .2-.3" = list(X = Xb, TR = .8, f1 = .2, f2 = .3),
    "Rs = 40" = list(X = Xb, TR = .72, f1 = .31, f2 = .43, Rs = 40),
    "single column" = list(X = Xb[, 1, drop = FALSE], TR = .72, f1 = .31, f2 = .43),
    "plain numeric vector" = list(X = v, TR = .72, f1 = .31, f2 = .43),
    "f2 above Nyquist (both error)" = list(X = Xb, TR = .72, f1 = .31, f2 = .70),
    "f1 >= f2 (both error)" = list(X = Xb, TR = .72, f1 = .43, f2 = .31)
  ), 1e-10)
  tsx <- stats::ts(rnorm(100))
  r_old <- run(OLD$bandstop_filter(tsx, .72, .31, .43)); r_new <- run(NEW$bandstop_filter(tsx, .72, .31, .43))
  expect_diff("bandstop_filter ts input: new runs and returns 100 rows",
    !r_new$error && NROW(r_new$value) == 100, paste("old:", if (r_old$error) errmsg(r_old) else "ran (length check was vacuous)"))
  if (!r_old$error && !r_new$error) {
    record("bandstop_filter ts input: values agree when old runs", "SAME",
           if (same(as.numeric(r_old$value), as.numeric(r_new$value), 1e-10)) "PASS" else "FAIL")
  }
} else skip("bandstop_filter", "gsignal not installed")


# =============================================================================
# Downstream: unchanged high-level function that calls changed code
# =============================================================================
if (is.function(OLD$norm_BOLD) && is.function(NEW$norm_BOLD)) {
  covers("norm_BOLD")
  set.seed(12)
  nV <- 120; nT <- 80
  Bn <- matrix(rnorm(nV * nT, 1000, 10), nV, nT)
  nuis <- cbind(sin(seq_len(nT) / 5))
  cmp_grid("norm_BOLD (downstream)", "norm_BOLD", list(
    "no filter" = list(BOLD = Bn, hpf = NULL, scale_by = "none"),
    "hpf" = list(BOLD = Bn, TR = 1, hpf = .01, scale_by = "none"),
    "nuisance" = list(BOLD = Bn, nuisance = nuis, hpf = NULL, scale_by = "none"),
    "drop_first" = list(BOLD = Bn, drop_first = 5, hpf = NULL, scale_by = "none"),
    "scale_by mean" = list(BOLD = Bn, hpf = NULL, scale_by = "mean", scale_sm_FWHM = 0),
    "scale_by sd" = list(BOLD = Bn, hpf = NULL, scale_by = "sd", scale_sm_FWHM = 0),
    "center_cols" = list(BOLD = Bn, hpf = NULL, center_cols = TRUE, scale_by = "none")
  ), 1e-8)
  a <- list(BOLD = Bn, TR = 1, hpf = .01, scale_by = "none")
  cmp_grid("norm_BOLD lpf only",   "norm_BOLD", list(lpf   = c(a, lpf = .1)), 1e-8)
  cmp_grid("norm_BOLD scrub only", "norm_BOLD", list(scrub = c(a, list(scrub = c(10, 30, 31)))), 1e-8)
  cmp_grid("norm_BOLD (downstream)", "norm_BOLD", list(
    "hpf + lpf + scrub" = list(BOLD = Bn, TR = 1, hpf = .01, lpf = .1, scrub = c(10, 30, 31), scale_by = "none")
  ), 1e-8)
} else skip("norm_BOLD (downstream)", "norm_BOLD not found in both OLD and NEW")


# =============================================================================
# Coverage: is everything that changed between OLD and NEW exercised above?
# =============================================================================
cat(sprintf("\nOLD: %s\nNEW: %s\n", normalizePath(old_dir), normalizePath(new_dir)))
cat(sprintf("Functions: %d modified, %d added, %d removed\n", length(fn_modified), length(fn_added), length(fn_removed)))
if (length(fn_modified)) cat("  modified:", paste(sort(fn_modified), collapse = ", "), "\n")
if (length(fn_added))    cat("  added:   ", paste(sort(fn_added), collapse = ", "), "\n")
if (length(fn_removed))  cat("  removed: ", paste(sort(fn_removed), collapse = ", "), "\n")

for (n_ in setdiff(c(fn_modified, fn_added), covered)) {
  record(n_, "COVERAGE", "GAP", if (n_ %in% fn_added) "added function; no check above" else "modified function; no check above")
}
# Unchanged functions that call changed ones (argument-validation helpers are
#  called everywhere; their changes are checked directly above).
changed_calls <- setdiff(c(fn_modified, fn_added), c("is_1", "is_posNum", "is_integer", "match_exactly"))
for (n_ in setdiff(fn_new, c(fn_modified, fn_added, covered))) {
  b <- tryCatch(all.names(body(get(n_, envir = NEW))), error = function(e) character())
  hit <- intersect(changed_calls, b)
  if (length(hit)) record(n_, "COVERAGE", "INFO", paste("unchanged; calls changed:", paste(hit, collapse = ", ")))
}


# =============================================================================
# Report
# =============================================================================
res <- do.call(rbind, results)
res$detail <- substr(res$detail, 1, 110)
cat("\n================ SUMMARY ================\n")
print(table(res$kind, res$status))
bad <- res[res$status == "FAIL", ]
if (nrow(bad) > 0) {
  cat("\n---------------- FAILURES ----------------\n")
  print(bad[, c("kind", "label", "detail")], row.names = FALSE, right = FALSE)
} else {
  cat("\nNo failures.\n")
}
cov_rows <- res[res$kind == "COVERAGE", ]
if (nrow(cov_rows) > 0) {
  cat("\n---------------- COVERAGE ----------------\n")
  print(cov_rows[, c("status", "label", "detail")], row.names = FALSE, right = FALSE)
}
if (any(res$status == "SKIP")) {
  cat("\n---------------- SKIPPED -----------------\n")
  print(res[res$status == "SKIP", c("label", "detail")], row.names = FALSE, right = FALSE)
}
out_csv <- file.path(script_dir, "compare_0.9.1_denoising_results.csv")
utils::write.csv(res, out_csv, row.names = FALSE)
cat("\nFull table written to", out_csv, "\n")
if (!interactive() && nrow(bad) > 0) quit(status = 1)
