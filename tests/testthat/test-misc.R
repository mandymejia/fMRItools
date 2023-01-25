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

})
