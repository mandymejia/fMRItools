check_wb <- function() {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }
}

test_that("Miscellaneous functions are working", {
  check_wb()

  tdir <- tempdir()

  # Do the tests
  mat <- matrix(rnorm(60), nrow=10)
  colCenter(mat)
  hat_matrix(mat)
})
