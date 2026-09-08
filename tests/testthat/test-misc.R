test_that("Miscellaneous functions are working", {

  tdir <- tempdir()

  # Do the tests

  # Other -----
  dmat <- cbind(matrix(rnorm(60), nrow=10), 1)
  dmat <- scale_design_mat(dmat)
  testthat::expect_equal(
    max(colMeans(dmat)), 0
  )
  dmat <- validate_design_mat(dmat)
  myhat <- hat_matrix(dmat)
  testthat::expect_equal(
    max(abs(myhat %*% myhat - myhat)), 0
  )

  set.seed(0)
  mat <- matrix(runif(100, min=-1, max=1), nrow=10)
  mat[upper.tri(mat)] <- mat[lower.tri(mat)]
  plot_FC(mat)
  plot_FC(mat, zlim=c(-.8, .8), diag_val=1, title="ABC", cleg_ticks_by=.2, lines="all", lines_col="white")
  plot_FC_gg(mat)
  plot_FC_gg(mat, title="ABC", legTitle="DEF", lim=c(-.8, .8), diagVal=0)
  plot_FC_gg(mat, labs=as.character(seq(10)))
  plot_FC_gg(mat, group_divs=c(1,2,5,8,10))
  plot_FC_gg(mat, group_divs=c(2,10), uppertri_means=FALSE, divColor="white", labs=as.character(seq(10)))
  plot_FC_gg(mat, group_divs=c(1,2,5,8,10), labs=c("a", "b", "c", "d", "e"))

  # `norm_BOLD` and `dual_reg`, `dual_reg_parc` -----
  nT <- 70
  nV <- 4000
  nQ <- 13
  mU <- matrix(rnorm(nV*nQ), nrow=nV)
  mS <- mU %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ)
  BOLD <- mS + rnorm(nV*nT, sd=.1) + 11

  testthat::expect_equal(
    BOLD,
    norm_BOLD(BOLD, center_rows=FALSE, scale_by="none", hpf=0)
  )

  nBOLD_mean <- norm_BOLD(BOLD, scale_sm_FWHM=0, hpf=NULL)
  testthat::expect_all_true(rowMeans(nBOLD_mean) < 1e-8)

  nBOLD_sd <- norm_BOLD(BOLD, scale_by="sd", scale_sm_FWHM=0, hpf=NULL)
  testthat::expect_all_equal(matrixStats::rowVars(nBOLD_sd), 1)

  testthat::expect_all_equal(abs(diag(
    cor(t(nBOLD_mean[seq(33,50),]), t(nBOLD_sd[seq(33,50),]))
  )), 1)

  if (is.null(ciftiTools:::ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }

  surfL <- ciftiTools::load_surf(resamp_res=3000)
  xii <- convert_to_dscalar(resample_xifti(load_parc(), resamp_res=3000)) + 5
  nV <- nrow(surfL$vertices)
  BOLD2 <- BOLD[seq(nV),] + c(as.matrix(xii$data$cortex_left))
  xiiL <- as.xifti(surfL=surfL)

  nBOLD_mean_0 <- norm_BOLD(BOLD2, hpf=NULL, scale_sm_FWHM=0)
  nBOLD_mean_sm <- norm_BOLD(BOLD2, hpf=NULL, scale_sm_xifti=xiiL)
  nBOLD_mean_Inf <- norm_BOLD(BOLD2, hpf=NULL, scale_sm_FWHM=Inf)

  bvars <- xiiL
  bvars$data$cortex_left <- cbind(
    rowVars(nBOLD_mean_0), rowVars(nBOLD_mean_sm), rowVars(nBOLD_mean_1)
  )
  plot(bvars, idx=seq(3))

  nBOLD_sd_0 <- norm_BOLD(BOLD2, scale_by="sd", hpf=NULL, scale_sm_FWHM=0)
  nBOLD_sd_sm <- norm_BOLD(BOLD2, scale_by="sd", hpf=NULL, scale_sm_xifti=xiiL)
  nBOLD_sd_Inf <- norm_BOLD(BOLD2, scale_by="sd", hpf=NULL, scale_sm_FWHM=Inf)

  bvars <- xiiL
  bvars$data$cortex_left <- cbind(
    rowVars(nBOLD_sd_0), rowVars(nBOLD_sd_sm), rowVars(nBOLD_sd_Inf)
  )
  plot(bvars, idx=seq(3))

})
