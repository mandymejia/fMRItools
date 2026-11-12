# Make internal (non-exported) functions available when this file is
# source()d from the runner script rather than run inside the package.
for (.nm in setdiff(ls(asNamespace("fMRItools"), all.names=TRUE), ls(globalenv()))) {
  assign(.nm, get(.nm, envir=asNamespace("fMRItools")), envir=globalenv())
}
rm(.nm)

# Build --> Install and Restart
# [Edit this] path to the Workbench for your computer.
# [TO DO]: shorten to just "~/Applications" after ciftiTools update
my_wb <- "~/Applications/ConnectomeWorkbench.app/Contents/usr/bin"

library(testthat)
library(fMRItools)
library(ciftiTools)

if (interactive()) { ciftiTools.setOption("wb_path", my_wb) }
tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }

source(file.path(tests_dir, "test-misc.R"))
source(file.path(tests_dir, "test-norm_BOLD.R"))
source(file.path(tests_dir, "test-varDecomp.R"))
source(file.path(tests_dir, "test-preprocessing.R"))
source(file.path(tests_dir, "test-CompCor.R"))
