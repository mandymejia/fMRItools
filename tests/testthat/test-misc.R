test_that("Miscellaneous functions are working", {
  tdir <- tempdir()

  # Do the tests
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
})
