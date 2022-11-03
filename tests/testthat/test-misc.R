test_that("Miscellaneous functions are working", {
  tdir <- tempdir()

  # Do the tests
  mat <- matrix(rnorm(60), nrow=10)
  hat_matrix(mat)
})
