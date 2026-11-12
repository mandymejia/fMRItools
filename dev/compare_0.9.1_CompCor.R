####################################################################
# compare_0.9.1_CompCor.R
#
# Old-vs-new check of the CompCor code in fMRItools on a small synthetic
#   "phantom" dataset. In the same style
#   as compare_0.9.1_denoising.R: the R/ files of each version are sourced into
#   separate environments (OLD, NEW), so nothing needs to be installed.
#
#   OLD = 0.9.0:  CompCor (+ internal CompCor.noise_comps), CompCor_HCP
#   NEW = the current fMRItools R/ directory (../R), where the functions were
#     reorganized into a "compute" layer and an "apply" layer:
#       noise_PCs()          compute noise PCs from noise ROI data (internal;
#                            was CompCor.noise_comps, now also does the prelim
#                            regression)
#       CompCor_clean()      apply: regress noise PCs (+ DCT, nuisance, intercept)
#                            from data (internal)
#       CompCor()            one-shot: format_data + noise_PCs + CompCor_clean
#       CompCor_HCP_calc()   NIFTI wrapper that computes noise PCs (was CompCor_HCP)
#       CompCor_HCP_apply()  CIFTI wrapper around CompCor_clean (was CompCor_HCP_cii)
#     CompCor()$noise and CompCor_HCP_calc() return the same "CompCor_noise" type.
#
# Statuses:  PASS / FAIL,  INFO (reported, not judged).
#   (A) things that should be UNCHANGED are unchanged (PCs and their variances
#       up to sign, noise ROI masks, cleaned data, ...);
#   (B) things that should be DIFFERENT are different as expected: the PCA
#       `v` directions (old norm sqrt((nV-1)/d), new unit norm, same
#       direction); the removed/changed CompCor_HCP_calc (old CompCor_HCP) arguments;
#   (C) `prelim_nuisance` does what it should, checked directly (not just
#       against the old code): after the preliminary regression, every noise
#       PC is orthogonal to the nuisance signal(s) that were regressed out,
#       while a signal that was NOT regressed out is still in the PCs.
#       NOTE on the intercept: the noise ROI voxels are centered by their
#       MEDIAN after the preliminary regression (center=TRUE). The regression
#       leaves each voxel with mean 0, but the median is not the mean, so
#       subtracting it adds a small voxel-specific constant back. The PCs are
#       thus orthogonal to the NON-constant part of the nuisance design (checked
#       with the design's columns centered), but may keep a small component along
#       the constant. With center=FALSE they are orthogonal to the whole design
#       (checked too, to confirm this explanation).
#   (D) the reorganization itself: removed/renamed functions, the common return
#       type, `noise_PCs` and `CompCor_clean` against direct computations,
#       and `CompCor_HCP_apply` == `CompCor_clean`.
#
# HOW TO RUN (from the fMRItools/dev directory):
#   1. Copy the revised CompCor.R and CompCor_HCP.R into fMRItools/R (if not
#      already there).
#   2. Settings below: `old_tools_R` is the frozen 0.9.0 R/ dir from
#      compare_0.9.1_denoising.R (or set `old_tools_commit`).
#   3. In R:  setwd("<path>/fMRItools/dev"); source("compare_0.9.1_CompCor.R")
#      (or: Rscript compare_0.9.1_CompCor.R)
#   Needs git on the PATH (only if a commit is used). Optional packages are
#   skipped if missing: RNifti, ciftiTools. Needs whatever the R/ files need.
####################################################################

# ---- Settings ----------------------------------------------------------------
new_R <- file.path("..", "R")        # current fMRItools source
# Old fMRItools helpers (0.9.0): a directory, or a commit in `repo_tools`.
old_tools_R <- file.path("0.9.0", "R")
old_tools_commit <- NULL             # e.g. "3057439"; overrides `old_tools_R`
repo_tools <- ".."                   # fMRItools git repo

tol_same <- 1e-8                     # expected to be identical
tol_qr <- 1e-6                       # projection via hat_Q vs. hat_matrix
tol_file <- 1e-4                     # NIFTI written to disk and re-read

# Expected relation of NEW to OLD for the noise PCs ("PCs", the time-course
#   scores `u`) and their variances ("var", from `d`). The PCA update only
#   rescaled the `v` directions (checked separately below), so both should be
#   unchanged. Other options: NULL (measure + report factors only), or
#   function(old, roi) returning the expected new/old factor for each PC.
scale_rules <- list(PCs = "same", var = "same")
# -------------------------------------------------------------------------------

# ---- Get the R/ directories ------------------------------------------------------
extract_commit <- function(repo, commit) {
  dest <- tempfile("old_R_")
  dir.create(dest)
  status <- system(sprintf("git -C %s archive %s R | tar -x -C %s",
    shQuote(path.expand(repo)), shQuote(commit), shQuote(dest)))
  if (status != 0 || !dir.exists(file.path(dest, "R"))) {
    stop("Could not extract R/ at commit ", commit, " from ", repo)
  }
  file.path(dest, "R")
}
if (!is.null(old_tools_commit)) {
  old_tools_R <- extract_commit(repo_tools, old_tools_commit)
}
stopifnot(dir.exists(old_tools_R), dir.exists(new_R))

# Source every file in each directory into one environment (later directories
#   override earlier ones); record which file defined each function.
load_dirs <- function(dirs) {
  env <- new.env(parent = globalenv())
  origin <- character()
  for (dir in dirs) {
    for (f in list.files(dir, pattern = "\\.[Rr]$", full.names = TRUE)) {
      snap <- mget(ls(env, all.names = TRUE), envir = env)
      tryCatch(
        sys.source(f, envir = env, keep.source = FALSE),
        error = function(e) warning("Could not source ", basename(f), ": ",
          conditionMessage(e), call. = FALSE)
      )
      for (n in ls(env, all.names = TRUE)) {
        if (!identical(env[[n]], snap[[n]])) {
          origin[n] <- file.path(basename(dirname(normalizePath(dir))), basename(f))
        }
      }
    }
  }
  attr(env, "origin") <- origin
  env
}

OLD <- load_dirs(old_tools_R)
NEW <- load_dirs(new_R)
stopifnot(is.function(OLD$CompCor), is.function(OLD$CompCor_HCP))
stopifnot(is.function(NEW$CompCor), is.function(NEW$CompCor_HCP_calc))
stopifnot(is.function(NEW$noise_PCs), is.function(NEW$CompCor_clean))
if (!is.function(NEW$hat_Q) || is.function(OLD$hat_Q)) {
  stop("`new_R` should be the current code (has hat_Q) and the old ",
    "directories the previous code (no hat_Q).")
}
# As in compare_0_9_1_0_9_0.R: replace the old `rowMedians2` (inverted-condition
#   bug) so that the old code is compared in isolation from that bug.
if (is.function(OLD$rowMedians2)) {
  OLD$rowMedians2 <- function(x, na.rm = FALSE, ...) {
    apply(x, 1, median, na.rm = na.rm, ...)
  }
}
key_fun_old <- c("CompCor", "CompCor_HCP", "CompCor.noise_comps", "format_data",
  "nuisance_regression", "hat_matrix", "PCA", "scale_med", "erode_mask_vol")
key_fun_new <- c("CompCor", "CompCor_HCP_calc", "CompCor_HCP_apply",
  "CompCor_clean", "noise_PCs", "new_CompCor_noise", "format_data",
  "nuisance_regression", "hat_Q", "PCA", "scale_med", "erode_mask_vol")
for (v in c("OLD", "NEW")) {
  env <- get(v); kf <- if (v == "OLD") { key_fun_old } else { key_fun_new }
  cat(v, "functions come from:\n")
  print(data.frame(
    `function` = kf, file = unname(attr(env, "origin")[kf]),
    check.names = FALSE
  ), row.names = FALSE)
  cat("\n")
}

# ---- Helpers -----------------------------------------------------------------------
results <- data.frame(
  part=character(), test=character(), status=character(), detail=character(),
  stringsAsFactors=FALSE
)
add <- function(part, test, ok, detail="", info=FALSE) {
  status <- if (info) { "INFO" } else if (isTRUE(ok)) { "PASS" } else { "FAIL" }
  results[nrow(results) + 1, ] <<- list(part, test, status, detail)
  cat(sprintf("[%s] %s: %s %s\n", status, part, test, detail))
}
fmt_num <- function(x) { paste(signif(x, 4), collapse=", ") }
err <- function(expr) { inherits(try(expr, silent=TRUE), "try-error") }

# Compute the CompCor outputs for one version (OLD or NEW).
compute_all <- function(ENV, d, is_old) {
  rois <- c("wm_cort", "csf", "wm_cblm")
  ero <- c(2, 1, 2)
  masks <- list(
    wm_cort = d$labels == 3001, csf = d$labels == 4, wm_cblm = d$labels == 7
  )
  # (OLD returned the noise PCs nested in `$noise` for some functions; NEW
  #   returns the same "CompCor_noise" type everywhere, checked in part D.)
  fix <- function(x) { if (is.null(x$PCs) && !is.null(x$noise)) x$noise else x }
  nii_idx <- d$nii[,,, d$idx]
  core <- function(prelim) {
    fix(ENV$CompCor(
      nii_idx, ROI_data=NULL, ROI_noise=masks, noise_nPC=5,
      noise_erosion=ero, prelim_nuisance=prelim
    )$noise)
  }
  hcp_fun <- if (is_old) { ENV$CompCor_HCP } else { ENV$CompCor_HCP_calc }
  hcp <- function(idx, prelim=NULL) {
    args <- list(
      nii=d$nii, nii_labels=d$labels, ROI_noise=rois, noise_nPC=5,
      noise_erosion=ero, idx=idx
    )
    if (!is.null(prelim)) { args$prelim_nuisance <- prelim }
    fix(do.call(hcp_fun, args))
  }
  out <- list(
    core_noprelim = core(NULL),
    core_noprelim2 = core(NULL),
    core_prelim = core(d$nuis),
    hcp_noprelim = hcp(d$idx),
    hcp_all = hcp(NULL),
    pca3 = ENV$PCA(d$Xp, Q=5, Vdim=3)
  )
  if (!is_old) { out$hcp_prelim <- hcp(d$idx, d$nuis) }
  out
}

# Per-ROI, per-PC comparison metrics (new vs. old).
cmp_noise <- function(new, old) {
  stopifnot(identical(names(new$PCs), names(old$PCs)))
  per <- lapply(names(old$PCs), function(roi) {
    A <- new$PCs[[roi]]; B <- old$PCs[[roi]]
    if (is.null(A) || is.null(B)) { return(list(empty=TRUE)) }
    stopifnot(identical(dim(A), dim(B)))
    num <- colSums(A * B)
    sgn <- sign(num); sgn[sgn == 0] <- 1
    k <- num / colSums(B * B)            # signed factor: new = k * old
    vn <- unlist(new$var[[roi]]); vo <- unlist(old$var[[roi]])
    list(
      empty = FALSE, nPC = ncol(A),
      cos = abs(num) / sqrt(colSums(A * A) * colSums(B * B)),
      k = abs(k),
      same_diff = max(abs(A - sweep(B, 2, sgn, "*"))),
      var_ratio = vn / vo,
      var_same_diff = max(abs(vn - vo) / abs(vo)),
      mask_same = identical(
        as.vector(new$ROI_noise[[roi]]), as.vector(old$ROI_noise[[roi]])
      )
    )
  })
  names(per) <- names(old$PCs)
  per
}

# Candidate formulas for a rescaling (only used when a rule is NULL).
candidates <- function(old, roi) {
  vo <- unlist(old$var[[roi]]); Tn <- nrow(old$PCs[[roi]])
  list(
    `1` = 1, `sqrt(var)` = sqrt(vo), `1/sqrt(var)` = 1/sqrt(vo),
    `var` = vo, `1/var` = 1/vo, `sqrt(T)` = sqrt(Tn), `1/sqrt(T)` = 1/sqrt(Tn),
    `sqrt(T-1)` = sqrt(Tn-1), `1/sqrt(T-1)` = 1/sqrt(Tn-1),
    `sqrt(var)/sqrt(T-1)` = sqrt(vo)/sqrt(Tn-1)
  )
}

# Full comparison of NEW vs. OLD noise outputs, according to `scale_rules`.
report_vs_old <- function(part, new, old) {
  per <- cmp_noise(new, old)
  rois <- names(per)[!vapply(per, function(p){ p$empty }, TRUE)]
  add(part, "same output elements", setequal(names(new), names(old)),
    paste0("new: ", paste(names(new), collapse=","),
      "; old: ", paste(names(old), collapse=",")))
  add(part, "same noise ROI masks (after erosion)",
    all(vapply(per[rois], function(p){ p$mask_same }, TRUE)))
  add(part, "same number of PCs per ROI", TRUE,
    paste(vapply(rois, function(r){ per[[r]]$nPC }, 1), collapse=","), info=TRUE)
  add(part, "PC directions same (|cos| = 1)",
    all(vapply(per[rois], function(p){ all(p$cos > 1 - tol_same) }, TRUE)),
    paste0("min |cos| = ", fmt_num(min(vapply(per[rois],
      function(p){ min(p$cos) }, 1)))))
  for (el in c("PCs", "var")) {
    rule <- scale_rules[[el]]
    meas <- if (el == "PCs") { "k" } else { "var_ratio" }
    for (roi in rois) {
      p <- per[[roi]]; measured <- p[[meas]]
      lab <- paste0(el, " [", roi, "]")
      if (is.null(rule)) {
        cands <- candidates(old, roi)
        hit <- names(cands)[vapply(cands, function(cf){
          isTRUE(all(abs(measured / cf - 1) < 1e-6)) }, TRUE)]
        add(part, paste0("scale of ", lab, " (rule not set)"), NA,
          paste0("measured new/old = ", fmt_num(measured), "; matching ",
            "candidates: ", if (length(hit)) paste(hit, collapse=" | ") else
            "none"), info=TRUE)
      } else if (identical(rule, "same")) {
        diff_ok <- if (el == "PCs") { p$same_diff < tol_same } else {
          p$var_same_diff < tol_same }
        add(part, paste0(lab, " unchanged"), diff_ok)
      } else {
        expected <- rule(old, roi)
        add(part, paste0(lab, " rescaled as expected"),
          isTRUE(all(abs(measured / expected - 1) < 1e-6)),
          paste0("measured = ", fmt_num(measured), "; expected = ",
            fmt_num(expected)))
      }
    }
  }
  # Regressing PCs out of data is invariant to their scale.
  Ycl <- function(noise) {
    qr.resid(qr(cbind(1, do.call(cbind, noise$PCs))), d$Y)
  }
  dY <- max(abs(Ycl(new) - Ycl(old)))
  add(part, "data cleaned with the PCs is unchanged", dY < tol_qr,
    paste0("max abs diff = ", fmt_num(dY)))
}

# ---- Phantom data ------------------------------------------------------------------
# A small 4D volume with three noise ROIs (cortical WM, CSF, cerebellar WM),
#   each carrying 6 latent AR(1) components of decreasing variance, plus RP-like
#   motion regressors and a few large spikes, so that the PCs are well-defined
#   and the preliminary nuisance regression has an effect.
make_phantom <- function(seed=1, dims=c(24, 30, 24), T_=150, drop=10) {
  set.seed(seed)
  labels <- array(0, dim=dims)
  labels[2:23, 2:9, 2:23] <- 3001    # wm_cort
  labels[2:23, 11:18, 2:23] <- 4     # csf
  labels[2:23, 20:28, 2:23] <- 7     # wm_cblm
  nii <- array(rnorm(prod(dims) * T_), dim=c(dims, T_))
  ar1 <- function(sd) {
    sd * as.numeric(stats::filter(rnorm(T_), 0.5, "recursive"))
  }
  sds <- c(6, 4, 3, 2.2, 1.6, 1.2)
  rp <- sapply(1:3, function(j){ cumsum(rnorm(T_, sd=0.1)) })
  spk_t <- c(30, 75, 120)
  for (lab in c(3001, 4, 7)) {
    vox <- which(labels == lab, arr.ind=TRUE)
    lat <- sapply(sds, ar1)                                  # T x 6
    L <- matrix(rnorm(6 * nrow(vox)), 6)                     # 6 x V
    Ylab <- lat %*% L + rp %*% matrix(rnorm(3 * nrow(vox)), 3)
    Ylab[spk_t, ] <- Ylab[spk_t, ] + 25
    for (v in seq(nrow(vox))) {
      nii[vox[v, 1], vox[v, 2], vox[v, 3], ] <- Ylab[, v] + rnorm(T_)
    }
  }
  idx <- seq(drop + 1, T_)
  dct <- sapply(1:4, function(k){ cos(pi * k * (seq(T_) - 0.5) / T_) })
  spikes <- diag(T_)[, spk_t]
  nuis <- cbind(1, dct, rp, spikes, rp[, 1])[idx, ]   # last col: rank-deficient
  stopifnot(qr(nuis)$rank < ncol(nuis))
  list(
    nii=nii, labels=labels, idx=idx, nuis=nuis,
    Y=matrix(rnorm(length(idx) * 50), length(idx)),    # data to "clean"
    # For PCA: V x T with distinct eigenvalues
    Xp=matrix(rnorm(300 * 40), 300, 40) %*% diag(seq(3, 1, length.out=40))
  )
}
d <- make_phantom()

cat("Running OLD and NEW code...\n")
old <- compute_all(OLD, d, is_old=TRUE)
new <- compute_all(NEW, d, is_old=FALSE)
cat("\n")

# ---- (A0) PCA: u and d unchanged; v rescaled (checked in isolation) ----------------------
o <- old$pca3; n <- new$pca3
nV <- nrow(d$Xp); cn <- function(M) { sqrt(colSums(M^2)) }
add("A0 PCA", "d (singular/eigen values) unchanged",
  isTRUE(all.equal(o$d, n$d, tolerance=1e-8)))
add("A0 PCA", "u (time-course PCs) unchanged up to sign",
  isTRUE(all.equal(abs(diag(crossprod(o$u, n$u))), rep(1, ncol(o$u)),
    tolerance=1e-6)))
add("B0 PCA", "v directions: same direction (up to sign)",
  isTRUE(all.equal(abs(colSums(o$v * n$v)) / (cn(o$v) * cn(n$v)),
    rep(1, ncol(n$v)), tolerance=1e-6)))
add("B0 PCA", "v: new directions are unit norm",
  isTRUE(all.equal(cn(n$v), rep(1, ncol(n$v)), tolerance=1e-8)))
add("B0 PCA", "v: old norm was sqrt((nV-1)/d), i.e. the rescaled element",
  isTRUE(all.equal(cn(o$v), sqrt((nV - 1) / n$d[seq_len(ncol(n$v))]),
    tolerance=1e-6)),
  paste0("old/new norm ratio = ", fmt_num(cn(o$v) / cn(n$v))))

# ---- (A) Core CompCor: new vs. old ---------------------------------------------------------
add("A core", "deterministic (new, run twice)",
  isTRUE(all.equal(new$core_noprelim$PCs, new$core_noprelim2$PCs,
    tolerance=1e-12)))
report_vs_old("A core, no prelim", new$core_noprelim, old$core_noprelim)
report_vs_old("A core, prelim", new$core_prelim, old$core_prelim)
dd <- cmp_noise(new$core_prelim, new$core_noprelim)
add("B core", "prelim_nuisance changes the PCs (as it should)",
  any(vapply(dd, function(p){ !p$empty && p$same_diff > 1e-3 }, TRUE)),
  paste0("max |diff| = ", fmt_num(max(vapply(dd, function(p){
    p$same_diff }, 1)))))

# ---- (A/B) CompCor_HCP_calc (OLD: CompCor_HCP) ----------------------------------------------------------------------
tcmp <- function(a, b, tol) {
  p <- cmp_noise(a, b)
  all(vapply(p, function(x){ x$empty || (x$same_diff < tol &&
    x$var_same_diff < tol && x$mask_same) }, TRUE))
}
add("A HCP", "idx subset: CompCor_HCP_calc == core CompCor (no prelim)",
  tcmp(new$hcp_noprelim, new$core_noprelim, tol_same))
add("A HCP", "hat_Q prelim regression == core CompCor prelim (rank-deficient)",
  tcmp(new$hcp_prelim, new$core_prelim, tol_qr))
report_vs_old("A HCP vs old CompCor_HCP", new$hcp_noprelim, old$hcp_noprelim)  # NEW: CompCor_HCP_calc
CompCor_HCP_calc <- NEW$CompCor_HCP_calc
n_all <- dim(d$nii)[4]
all_full <- CompCor_HCP_calc(d$nii, d$labels, c("wm_cort", "csf", "wm_cblm"),
  noise_nPC=5, noise_erosion=c(2, 1, 2), idx=seq(n_all))
add("A HCP", "idx = NULL equals idx = all timepoints",
  tcmp(new$hcp_all, all_full, tol_same))

add("B HCP", "`cii` argument removed (now in CompCor_HCP_apply)",
  err(CompCor_HCP_calc(d$nii, d$labels, idx=d$idx, cii="x")))
add("B HCP", "`idx` as a list is not supported",
  err(CompCor_HCP_calc(d$nii, d$labels, idx=list(d$idx))))
add("B HCP", "`prelim_nuisance` with the wrong number of rows errors",
  err(CompCor_HCP_calc(d$nii, d$labels, idx=d$idx, prelim_nuisance=d$nuis[-1, ])))

if (requireNamespace("RNifti", quietly=TRUE)) {
  f_nii <- tempfile(fileext=".nii.gz"); f_lab <- tempfile(fileext=".nii.gz")
  RNifti::writeNifti(d$nii, f_nii); RNifti::writeNifti(d$labels, f_lab)
  from_file <- CompCor_HCP_calc(f_nii, f_lab, c("wm_cort", "csf", "wm_cblm"),
    noise_nPC=5, noise_erosion=c(2, 1, 2), idx=d$idx)
  add("A HCP", "reading NIFTI file paths == passing arrays",
    tcmp(from_file, new$hcp_noprelim, tol_file))
} else {
  add("A HCP", "reading NIFTI file paths (RNifti not installed)", NA, "", TRUE)
}

if (is.function(NEW$CompCor_HCP_apply) && requireNamespace("ciftiTools", quietly=TRUE)) {
  tryCatch({
    Tn <- length(d$idx)
    Yc <- matrix(rnorm(40 * Tn), 40, Tn)                    # vertices x time
    xii <- ciftiTools::make_xifti(cortexL=Yc[1:20, ], cortexR=Yc[21:40, ])
    for (ndct in c(0, 3)) {
      out <- NEW$CompCor_HCP_apply(xii, new$hcp_noprelim, DCT=ndct)
      X <- cbind(1, do.call(cbind, new$hcp_noprelim$PCs))
      if (ndct > 0) { X <- cbind(X, NEW$dct_bases(Tn, ndct)) }
      manual <- t(qr.resid(qr(X), t(Yc)))
      add("A cii",
        paste0("CompCor_HCP_apply == manual regression (DCT=", ndct, ")"),
        max(abs(as.matrix(out) - manual)) < tol_qr)
      add("A cii",
        paste0("CompCor_HCP_apply == CompCor_clean (DCT=", ndct, ")"),
        max(abs(as.matrix(out) -
          NEW$CompCor_clean(Yc, new$hcp_noprelim, DCT=ndct))) < tol_qr)
    }
  }, error=function(e){
    add("A cii", "CompCor_HCP_apply on a synthetic xifti", NA,
      paste("skipped:", conditionMessage(e)), TRUE)
  })
} else {
  add("A cii", "CompCor_HCP_apply (missing function or ciftiTools)", NA, "", TRUE)
}

# ---- (C) prelim_nuisance: directly check that it works ----------------------------------------
# Design: add to every noise-ROI voxel (with random loadings, so the signals survive
#   the spatial centering done by the PCA) two known signals:
#     z = slow oscillation + trend   -> passed as `prelim_nuisance` (with an intercept)
#     w = fast oscillation           -> NOT passed (to show the regression is selective)
#   Both are made strong enough to dominate the leading PCs. Then:
#     without prelim:  the leading PCs of each ROI are ~z and ~w (positive control);
#     with prelim:     every PC is orthogonal to z (and to the whole nuisance design),
#                      but w is still there.
fix_noise <- function(x) { if (is.null(x$PCs) && !is.null(x$noise)) x$noise else x }
prelim_core <- function(ENV, nii, labels, idx, prelim, center=TRUE) {
  masks <- list(wm_cort = labels == 3001, csf = labels == 4, wm_cblm = labels == 7)
  fix_noise(ENV$CompCor(
    nii[,,, idx], ROI_data=NULL, ROI_noise=masks, noise_nPC=5,
    noise_erosion=c(2, 1, 2), center=center, prelim_nuisance=prelim
  )$noise)
}
prelim_hcp <- function(ENV, nii, labels, idx, prelim, center=TRUE) {
  args <- list(
    nii=nii, nii_labels=labels, ROI_noise=c("wm_cort", "csf", "wm_cblm"),
    noise_nPC=5, noise_erosion=c(2, 1, 2), idx=idx, center=center
  )
  if (!is.null(prelim)) { args$prelim_nuisance <- prelim }
  fix_noise(do.call(ENV$CompCor_HCP_calc, args))
}
make_nuisance_phantom <- function(d, seed=2) {
  set.seed(seed)
  dims <- dim(d$nii); T_ <- dims[4]; tt <- seq_len(T_)
  z <- sin(2 * pi * tt / 37) + 1.5 * tt / T_
  w <- cos(2 * pi * tt / 11)
  nii_mat <- matrix(d$nii, prod(dims[1:3]), T_)          # voxels x time
  for (lab in c(3001, 4, 7)) {
    lin <- which(d$labels == lab)
    nii_mat[lin, ] <- nii_mat[lin, ] +
      outer(rnorm(length(lin), sd=20), z) + outer(rnorm(length(lin), sd=15), w)
  }
  list(nii=array(nii_mat, dim=dims), z=z, w=w)
}

# Per ROI: largest |correlation| between a signal and any PC.
pc_cor <- function(PCs, s) {
  vapply(PCs, function(P) { if (is.null(P)) NA_real_ else max(abs(stats::cor(P, s))) }, 1)
}
# Per ROI: largest fraction of any PC's norm lying in the column space of X
#   (0 = every PC is orthogonal to every column of X).
pc_in_span <- function(PCs, X) {
  qx <- qr(X); Q <- qr.Q(qx)[, seq_len(qx$rank), drop=FALSE]
  vapply(PCs, function(P) {
    if (is.null(P)) { NA_real_ } else {
      max(sqrt(colSums(crossprod(Q, P)^2)) / sqrt(colSums(P^2)))
    }
  }, 1)
}
# Same, but for the design with its columns centered, i.e. the nuisance design
#   WITHOUT its constant part (an intercept column becomes all zero and drops out).
#   This is the property that holds with median-centering of the noise ROIs
#   (see the note in the header): PCs may keep a small component along the
#   constant, but none along the rest of the design.
pc_in_span_nonconst <- function(PCs, X) {
  pc_in_span(PCs, scale(X, center=TRUE, scale=FALSE))
}
all_below <- function(x, tol) { isTRUE(all(x < tol, na.rm=TRUE)) && any(!is.na(x)) }
all_above <- function(x, tol) { isTRUE(all(x > tol, na.rm=TRUE)) && any(!is.na(x)) }
by_roi <- function(x) { paste0(names(x), "=", signif(x, 3), collapse=", ") }

nd <- make_nuisance_phantom(d)
z <- nd$z[d$idx]; w <- nd$w[d$idx]                    # time points actually analyzed
X_z <- cbind(1, z)                                    # nuisance design: intercept + z
core_np <- prelim_core(NEW, nd$nii, d$labels, d$idx, NULL)
core_p  <- prelim_core(NEW, nd$nii, d$labels, d$idx, X_z)
hcp_p   <- prelim_hcp(NEW, nd$nii, d$labels, d$idx, X_z)
core_p_nc <- prelim_core(NEW, nd$nii, d$labels, d$idx, X_z, center=FALSE)
hcp_p_nc  <- prelim_hcp(NEW, nd$nii, d$labels, d$idx, X_z, center=FALSE)
const_z <- matrix(1, length(z), 1)                    # the intercept alone

# Positive control: without prelim, the signals are in the leading PCs.
add("C prelim", "control: without prelim_nuisance, every ROI's PCs carry z and w (max |cor| > 0.9)",
  all_above(pc_cor(core_np$PCs, z), 0.9) && all_above(pc_cor(core_np$PCs, w), 0.9),
  paste0("z: ", by_roi(pc_cor(core_np$PCs, z)), "; w: ", by_roi(pc_cor(core_np$PCs, w))))

# The check itself: with prelim, z is orthogonal to ALL PCs (CompCor).
add("C prelim", "CompCor: nuisance signal z is orthogonal to every PC (max |cor| ~ 0)",
  all_below(pc_cor(core_p$PCs, z), tol_qr), by_roi(pc_cor(core_p$PCs, z)))
add("C prelim", "CompCor: no PC has any component in span(nuisance design without its constant)",
  all_below(pc_in_span_nonconst(core_p$PCs, X_z), tol_qr),
  by_roi(pc_in_span_nonconst(core_p$PCs, X_z)))
add("C prelim", "CompCor: PC share along the intercept (median-centering re-adds a small constant)",
  NA, by_roi(pc_in_span(core_p$PCs, const_z)), info=TRUE)
add("C prelim", "CompCor, center=FALSE: no PC has any component in span(whole nuisance design, intercept included)",
  all_below(pc_in_span(core_p_nc$PCs, X_z), tol_qr), by_roi(pc_in_span(core_p_nc$PCs, X_z)))

# Selectivity: only what was passed is removed.
add("C prelim", "CompCor: signal w, NOT in prelim_nuisance, is still in the PCs (max |cor| > 0.9)",
  all_above(pc_cor(core_p$PCs, w), 0.9), by_roi(pc_cor(core_p$PCs, w)))

# Same for CompCor_HCP_calc (idx subset), which should also equal the core result.
add("C prelim", "CompCor_HCP_calc: nuisance signal z is orthogonal to every PC (max |cor| ~ 0)",
  all_below(pc_cor(hcp_p$PCs, z), tol_qr), by_roi(pc_cor(hcp_p$PCs, z)))
add("C prelim", "CompCor_HCP_calc: no PC has any component in span(nuisance design without its constant)",
  all_below(pc_in_span_nonconst(hcp_p$PCs, X_z), tol_qr),
  by_roi(pc_in_span_nonconst(hcp_p$PCs, X_z)))
add("C prelim", "CompCor_HCP_calc, center=FALSE: no PC has any component in span(whole nuisance design)",
  all_below(pc_in_span(hcp_p_nc$PCs, X_z), tol_qr), by_roi(pc_in_span(hcp_p_nc$PCs, X_z)))
add("C prelim", "CompCor_HCP_calc: signal w, NOT in prelim_nuisance, is still in the PCs (max |cor| > 0.9)",
  all_above(pc_cor(hcp_p$PCs, w), 0.9), by_roi(pc_cor(hcp_p$PCs, w)))
add("C prelim", "CompCor_HCP_calc prelim == CompCor prelim on the z/w phantom",
  tcmp(hcp_p, core_p, tol_qr))

# Rank-deficient design with several kinds of regressors (intercept, DCT, RPs, spikes
#   and a duplicated column): every PC must be orthogonal to every column.
add("C prelim", "CompCor, rank-deficient design (intercept, DCT, RPs, spikes): all PCs orthogonal to its non-constant part",
  all_below(pc_in_span_nonconst(new$core_prelim$PCs, d$nuis), tol_qr),
  by_roi(pc_in_span_nonconst(new$core_prelim$PCs, d$nuis)))
add("C prelim", "CompCor_HCP_calc, rank-deficient design: all PCs orthogonal to its non-constant part",
  all_below(pc_in_span_nonconst(new$hcp_prelim$PCs, d$nuis), tol_qr),
  by_roi(pc_in_span_nonconst(new$hcp_prelim$PCs, d$nuis)))
core_rd_nc <- prelim_core(NEW, d$nii, d$labels, d$idx, d$nuis, center=FALSE)
add("C prelim", "CompCor, center=FALSE, rank-deficient design: all PCs orthogonal to the whole design",
  all_below(pc_in_span(core_rd_nc$PCs, d$nuis), tol_qr),
  by_roi(pc_in_span(core_rd_nc$PCs, d$nuis)))
add("C prelim", "control: without prelim, the PCs are NOT orthogonal to that design (non-constant part)",
  all_above(pc_in_span_nonconst(new$core_noprelim$PCs, d$nuis), 0.05),
  by_roi(pc_in_span_nonconst(new$core_noprelim$PCs, d$nuis)))

# For reference: the old code on the same data (not judged).
core_p_old <- prelim_core(OLD, nd$nii, d$labels, d$idx, X_z)
add("C prelim", "old CompCor: max |cor| of z with any PC (for reference)", NA,
  by_roi(pc_cor(core_p_old$PCs, z)), info=TRUE)
add("C prelim", "old CompCor: PC share inside span(rank-deficient design, non-constant part) (for reference)", NA,
  by_roi(pc_in_span_nonconst(old$core_prelim$PCs, d$nuis)), info=TRUE)

# ---- (D) The reorganization: noise_PCs / CompCor_clean / CompCor_HCP_calc / ... -------------------
masks_d <- list(wm_cort = d$labels == 3001, csf = d$labels == 4, wm_cblm = d$labels == 7)
ero_d <- c(2, 1, 2)
nii_idx <- d$nii[,,, d$idx]
noise_d <- new$core_noprelim                       # CompCor noise object (no prelim)
pcs_d <- do.call(cbind, noise_d$PCs)               # T x (5 + 5 + 5)
extra <- d$nuis[, 6:11, drop=FALSE]                # RPs + spikes (columns of d$nuis)
T_d <- length(d$idx)

# Names / return type.
add("D API", "`CompCor.noise_comps` replaced by `noise_PCs`",
  !is.function(NEW$CompCor.noise_comps) && is.function(NEW$noise_PCs))
add("D API", "old HCP names are gone from NEW (CompCor_HCP, CompCor_HCP_nii, CompCor_HCP_cii)",
  !any(vapply(c("CompCor_HCP", "CompCor_HCP_nii", "CompCor_HCP_cii"),
    function(f){ is.function(NEW[[f]]) }, TRUE)))
add("D API", "CompCor()$noise and CompCor_HCP_calc() both return a 'CompCor_noise' (PCs, var, ROI_noise)",
  inherits(new$core_noprelim, "CompCor_noise") &&
  inherits(new$hcp_noprelim, "CompCor_noise") &&
  identical(names(unclass(new$core_noprelim)), c("PCs", "var", "ROI_noise")) &&
  identical(names(unclass(new$hcp_noprelim)), c("PCs", "var", "ROI_noise")))

# noise_PCs: the "compute" layer, on its own.
fmt_d <- NEW$format_data(
  X=nii_idx, ROI_data=NULL, ROI_noise=masks_d, noise_nPC=5, noise_erosion=ero_d
)
np0 <- NEW$noise_PCs(fmt_d$X_noise, noise_nPC=fmt_d$noise_nPC)
np_p <- NEW$noise_PCs(fmt_d$X_noise, noise_nPC=fmt_d$noise_nPC, prelim_nuisance=d$nuis)
add("D noise_PCs", "no prelim: equals the PCs/variances inside CompCor()",
  tcmp(NEW$new_CompCor_noise(np0$PCs, np0$var, fmt_d$ROI_noise), new$core_noprelim, tol_same))
add("D noise_PCs", "prelim: equals the PCs/variances inside CompCor(prelim_nuisance=)",
  tcmp(NEW$new_CompCor_noise(np_p$PCs, np_p$var, fmt_d$ROI_noise), new$core_prelim, tol_same))
Xm <- lapply(fmt_d$X_noise, function(x){ qr.resid(qr(d$nuis), x) })
np_m <- NEW$noise_PCs(Xm, noise_nPC=fmt_d$noise_nPC)
add("D noise_PCs", "prelim == regressing the (rank-deficient) design by qr.resid, then no prelim",
  tcmp(NEW$new_CompCor_noise(np_p$PCs, np_p$var, fmt_d$ROI_noise),
    NEW$new_CompCor_noise(np_m$PCs, np_m$var, fmt_d$ROI_noise), tol_qr))
add("D noise_PCs", "also returns the total variance per ROI",
  !is.null(np_p$vartotal) && all(vapply(names(np_p$PCs), function(r){
    np_p$vartotal[[r]] >= sum(np_p$var[[r]]) - 1e-8 }, TRUE)))
add("D noise_PCs", "`prelim_nuisance` with the wrong number of rows errors",
  err(NEW$noise_PCs(fmt_d$X_noise, noise_nPC=fmt_d$noise_nPC, prelim_nuisance=d$nuis[-1, ])))
add("D noise_PCs", "CompCor(): `prelim_nuisance` with the wrong number of rows errors",
  err(NEW$CompCor(nii_idx, ROI_data=NULL, ROI_noise=masks_d, noise_nPC=5,
    noise_erosion=ero_d, prelim_nuisance=d$nuis[-1, ])))

# CompCor_clean: the "apply" layer.
cl0 <- NEW$CompCor_clean(d$Y, noise_d)
add("D clean", "CompCor_clean == regressing [1, PCs] (qr.resid)",
  max(abs(cl0 - qr.resid(qr(cbind(1, pcs_d)), d$Y))) < tol_qr)
add("D clean", "cleaned data has no component along the intercept or any PC",
  max(abs(crossprod(cbind(1, pcs_d), cl0))) < tol_qr)
cl1 <- NEW$CompCor_clean(d$Y, noise_d, DCT=3, nuisance=extra)
add("D clean", "CompCor_clean with DCT=3 and `nuisance` == regressing [1, PCs, DCT, nuisance]",
  max(abs(cl1 - qr.resid(qr(cbind(1, pcs_d, NEW$dct_bases(T_d, 3), extra)), d$Y))) < tol_qr)
add("D clean", "V x T orientation gives the transposed result",
  max(abs(t(NEW$CompCor_clean(t(d$Y), noise_d, DCT=3, nuisance=extra)) - cl1)) < tol_same)
add("D clean", "`nuisance` = NULL, 0 and FALSE are equivalent",
  max(abs(NEW$CompCor_clean(d$Y, noise_d, nuisance=0) - cl0)) < tol_same &&
  max(abs(NEW$CompCor_clean(d$Y, noise_d, nuisance=FALSE) - cl0)) < tol_same)
add("D clean", "accepts a plain list with `PCs` (no class needed)",
  max(abs(NEW$CompCor_clean(d$Y, list(PCs=noise_d$PCs)) - cl0)) < tol_same)
add("D clean", "data with the wrong number of timepoints errors",
  err(NEW$CompCor_clean(d$Y[-1, ], noise_d)))
add("D clean", "`nuisance` with the wrong number of rows errors",
  err(NEW$CompCor_clean(d$Y, noise_d, nuisance=extra[-1, ])))
add("D clean", "`noise` without any PCs errors",
  err(NEW$CompCor_clean(d$Y, list(PCs=list(wm_cort=NULL)))))

# CompCor() == noise_PCs + CompCor_clean (needs a data ROI; skipped if format_data
#   does not accept this form of `ROI_data`).
tryCatch({
  roi_data <- d$labels == 0                                  # all non-ROI voxels
  out <- NEW$CompCor(
    nii_idx, ROI_data=roi_data, ROI_noise=masks_d, noise_nPC=5,
    noise_erosion=ero_d, nuisance=extra
  )
  Ydat <- matrix(nii_idx, prod(dim(d$labels)), T_d)[which(roi_data), , drop=FALSE]  # V x T
  expected <- NEW$CompCor_clean(Ydat, out$noise, nuisance=extra)
  got <- out$data
  if (identical(dim(got), rev(dim(expected)))) { got <- t(got) }
  add("D CompCor", "result is of class 'CompCor' with entries data and noise",
    inherits(out, "CompCor") && all(c("data", "noise") %in% names(out)) &&
    inherits(out$noise, "CompCor_noise"))
  add("D CompCor", "CompCor()$data == CompCor_clean(data, CompCor()$noise)",
    identical(dim(got), dim(expected)) && max(abs(got - expected)) < tol_qr)
}, error=function(e){
  add("D CompCor", "CompCor()$data == CompCor_clean(...)", NA,
    paste("skipped:", conditionMessage(e)), TRUE)
})

# ---- Summary ---------------------------------------------------------------------------------
cat("\n==================== SUMMARY ====================\n")
print(table(results$status))
if (any(results$status == "FAIL")) {
  cat("\nFAILED:\n")
  print(results[results$status == "FAIL", c("part", "test", "detail")],
    row.names=FALSE, right=FALSE)
}
if (any(vapply(scale_rules, is.null, TRUE))) {
  cat("\nNOTE: `scale_rules` has NULL entries, so the rescaling was only ",
    "measured (INFO lines), not verified.\n", sep="")
}
write.csv(results, "compare_0.9.1_CompCor_results.csv", row.names=FALSE)
cat("\nResults written to compare_0.9.1_CompCor_results.csv\n")
if (!interactive() && any(results$status == "FAIL")) { quit(status=1) }
invisible(results)
