# Build --> Install and Restart

library(testthat)
library(fMRItools)

tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }

source(file.path(tests_dir, "test-misc.R"))
