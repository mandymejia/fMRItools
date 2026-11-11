####################################################################
# compare_0.9.1_CompCor.R
#
# Old-vs-new check of the CompCor code in fMRItools (CompCor, CompCor_HCP,
#   CompCor_HCP_cii), on a small synthetic "phantom" dataset. In the same style
#   as compare_0.9.1_denoising.R: the R/ files of each version are sourced into
#   separate environments (OLD, NEW), so nothing needs to be installed.
#
#   OLD = 0.9.0
#   NEW = the current fMRItools R/ directory (../R).
#
# Statuses:  PASS / FAIL,  INFO (reported, not judged).
#   (A) things that should be UNCHANGED are unchanged (PCs and their variances
#       up to sign, noise ROI masks, cleaned data, ...);
#   (B) things that should be DIFFERENT are different as expected: the PCA
#       `v` directions (old norm sqrt((nV-1)/d), new unit norm, same
#       direction); the removed/changed CompCor_HCP arguments.
#
# HOW TO RUN (from the fMRItools/dev directory):
#   1. Copy the revised CompCor_HCP.R into fMRItools/R (if not already there).
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
stopifnot(is.function(NEW$CompCor), is.function(NEW$CompCor_HCP))
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
key_fun <- c("CompCor", "CompCor_HCP", "CompCor.noise_comps", "format_data",
  "nuisance_regression", "hat_matrix", "PCA", "scale_med", "erode_mask_vol")
cat("OLD functions come from:\n")
print(data.frame(
  `function` = key_fun, file = unname(attr(OLD, "origin")[key_fun]),
  check.names = FALSE
), row.names = FALSE)
cat("\n")

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

# Compute the CompCor outputs for one version (OLD or NEW).
compute_all <- function(ENV, d, is_old) {
  rois <- c("wm_cort", "csf", "wm_cblm")
  ero <- c(2, 1, 2)
  masks <- list(
    wm_cort = d$labels == 3001, csf = d$labels == 4, wm_cblm = d$labels == 7
  )
  fix <- function(x) { if (is.null(x$PCs) && !is.null(x$noise)) x$noise else x }
  nii_idx <- d$nii[,,, d$idx]
  core <- function(prelim) {
    fix(ENV$CompCor(
      nii_idx, ROI_data=NULL, ROI_noise=masks, noise_nPC=5,
      noise_erosion=ero, prelim_nuisance=prelim
    )$noise)
  }
  hcp <- function(idx, prelim=NULL) {
    args <- list(
      nii=d$nii, nii_labels=d$labels, ROI_noise=rois, noise_nPC=5,
      noise_erosion=ero, idx=idx
    )
    if (!is.null(prelim)) { args$prelim_nuisance <- prelim }
    fix(do.call(ENV$CompCor_HCP, args))
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

# ---- (A/B) CompCor_HCP ----------------------------------------------------------------------
tcmp <- function(a, b, tol) {
  p <- cmp_noise(a, b)
  all(vapply(p, function(x){ x$empty || (x$same_diff < tol &&
    x$var_same_diff < tol && x$mask_same) }, TRUE))
}
add("A HCP", "idx subset: CompCor_HCP == core CompCor (no prelim)",
  tcmp(new$hcp_noprelim, new$core_noprelim, tol_same))
add("A HCP", "hat_Q prelim regression == core CompCor prelim (rank-deficient)",
  tcmp(new$hcp_prelim, new$core_prelim, tol_qr))
report_vs_old("A HCP vs old CompCor_HCP", new$hcp_noprelim, old$hcp_noprelim)
CompCor_HCP <- NEW$CompCor_HCP
n_all <- dim(d$nii)[4]
all_full <- CompCor_HCP(d$nii, d$labels, c("wm_cort", "csf", "wm_cblm"),
  noise_nPC=5, noise_erosion=c(2, 1, 2), idx=seq(n_all))
add("A HCP", "idx = NULL equals idx = all timepoints",
  tcmp(new$hcp_all, all_full, tol_same))

err <- function(expr) { inherits(try(expr, silent=TRUE), "try-error") }
add("B HCP", "`cii` argument removed (now in CompCor_HCP_cii)",
  err(CompCor_HCP(d$nii, d$labels, idx=d$idx, cii="x")))
add("B HCP", "`idx` as a list is not supported",
  err(CompCor_HCP(d$nii, d$labels, idx=list(d$idx))))
add("B HCP", "`prelim_nuisance` with the wrong number of rows errors",
  err(CompCor_HCP(d$nii, d$labels, idx=d$idx, prelim_nuisance=d$nuis[-1, ])))

if (requireNamespace("RNifti", quietly=TRUE)) {
  f_nii <- tempfile(fileext=".nii.gz"); f_lab <- tempfile(fileext=".nii.gz")
  RNifti::writeNifti(d$nii, f_nii); RNifti::writeNifti(d$labels, f_lab)
  from_file <- CompCor_HCP(f_nii, f_lab, c("wm_cort", "csf", "wm_cblm"),
    noise_nPC=5, noise_erosion=c(2, 1, 2), idx=d$idx)
  add("A HCP", "reading NIFTI file paths == passing arrays",
    tcmp(from_file, new$hcp_noprelim, tol_file))
} else {
  add("A HCP", "reading NIFTI file paths (RNifti not installed)", NA, "", TRUE)
}

if (is.function(NEW$CompCor_HCP_cii) && requireNamespace("ciftiTools", quietly=TRUE)) {
  tryCatch({
    Tn <- length(d$idx)
    Yc <- matrix(rnorm(40 * Tn), 40, Tn)                    # vertices x time
    xii <- ciftiTools::make_xifti(cortexL=Yc[1:20, ], cortexR=Yc[21:40, ])
    for (ndct in c(0, 3)) {
      out <- NEW$CompCor_HCP_cii(xii, new$hcp_noprelim, DCT=ndct)
      X <- cbind(1, do.call(cbind, new$hcp_noprelim$PCs))
      if (ndct > 0) { X <- cbind(X, NEW$dct_bases(Tn, ndct)) }
      manual <- t(qr.resid(qr(X), t(Yc)))
      add("A cii",
        paste0("CompCor_HCP_cii == manual regression (DCT=", ndct, ")"),
        max(abs(as.matrix(out) - manual)) < tol_qr)
    }
  }, error=function(e){
    add("A cii", "CompCor_HCP_cii on a synthetic xifti", NA,
      paste("skipped:", conditionMessage(e)), TRUE)
  })
} else {
  add("A cii", "CompCor_HCP_cii (missing function or ciftiTools)", NA, "", TRUE)
}

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
write.csv(results, "check_CompCor_pkg_results.csv", row.names=FALSE)
cat("\nResults written to check_CompCor_pkg_results.csv\n")
if (!interactive() && any(results$status == "FAIL")) { quit(status=1) }
invisible(results)
