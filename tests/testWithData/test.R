# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
# [Edit this] path to the Workbench.
ciftiTools.setOption("wb_path", "~/Desktop/workbench")

library(gifti)
library(RNifti)
library(oro.nifti)

# library(fMRItools)
roxygen2::roxygenize()

# file paths
# [Edit this] path to the data directory.
my_data <- "/Volumes/GoogleDrive/My Drive/MEJIA_LAB-Damon/Data"

subjects <- c(100307, 100408, 100610)
rs_cii_fnames <- file.path(my_data, "HCP-RestingState", c(
  paste0(subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
))
rs_nii_fnames <- file.path(my_data, "HCP-RestingState", c(
  paste0(subjects, "_rfMRI_REST1_LR.nii.gz"),
  paste0(subjects, "_rfMRI_REST2_LR.nii.gz")
))
rm(subjects)
rs_giiL_fnames <- list.files(
  file.path(my_data, "MSC-OneSubject"),
  "fsLR.sep.L.func.gii", full.names=TRUE
)
rs_giiR_fnames <- list.files(
  file.path(my_data, "MSC-OneSubject"),
  "fsLR.sep.R.func.gii", full.names=TRUE
)
sc_cii_fnames <- list.files(
  file.path(my_data, "HCP-Structural-dscalars"),
  full.names=TRUE
)
lb_gii_fnames <- list.files(file.path(my_data, "Label"), full.names=TRUE)
sf_gii_fnames <- list.files(
  file.path(my_data, "MSC-OneSubject/Surfs"), full.names=TRUE
)

# infer_format_ifti and as.matrix_ifti -----------------------------------------

# CIFTI
testthat::expect_equal(
  infer_format_ifti(rs_cii_fnames[1]),
  c("CIFTI", "dtseries")
)
z <- as.matrix_ifti(rs_cii_fnames[1], idx=1, meta=TRUE)
testthat::expect_equal(
  infer_format_ifti(read_xifti(rs_cii_fnames[1], idx=1)),
  c("xifti", "dtseries")
)
testthat::expect_equal(
  infer_format_ifti(ciftiTools.files()$cifti["dscalar"]),
  c("CIFTI", "dscalar")
)
testthat::expect_equal(
  infer_format_ifti(read_xifti(ciftiTools.files()$cifti["dscalar"])),
  c("xifti", "dscalar")
)
testthat::expect_equal(
  infer_format_ifti(ciftiTools.files()$cifti["dlabel"]),
  c("CIFTI", "dlabel")
)
testthat::expect_equal(
  infer_format_ifti(read_xifti(ciftiTools.files()$cifti["dlabel"])),
  c("xifti", "dlabel")
)
z <- as.matrix_ifti(read_cifti(ciftiTools.files()$cifti["dlabel"]), meta=TRUE)
testthat::expect_equal(
  infer_format_ifti(ciftiTools::template_xifti()),
  c("xifti", NA)
)

# GIFTI
testthat::expect_equal(
  infer_format_ifti(rs_giiL_fnames[1]),
  c("GIFTI", "metric")
)
z <- as.matrix_ifti(readgii(rs_giiR_fnames[1]), meta=TRUE)
testthat::expect_equal(
  infer_format_ifti(readgii(rs_giiL_fnames[1])),
  c("gifti", "metric")
)
testthat::expect_equal(
  infer_format_ifti(sf_gii_fnames[1]),
  c("GIFTI", "surf")
)
testthat::expect_error(as.matrix_ifti(sf_gii_fnames[1], meta=TRUE))

# CIFTI and GIFTI return the same result
z <- do.call(
  cbind,
  list(
    as.matrix_ifti(rs_giiL_fnames[1]),
    as.matrix_ifti(rs_giiR_fnames[1])
  )
)
z2 <- as.matrix_ifti(
  move_from_mwall(read_cifti(
    gsub("sep.L.func.gii", "dtseries.nii", rs_giiL_fnames[1])
  ), 0)
)
testthat::expect_equal(z, z2)
rm(z2)

# NIFTI
testthat::expect_equal(
  infer_format_ifti(rs_nii_fnames[1]),
  c("NIFTI", NA)
)
z <- as.matrix_ifti(rs_nii_fnames[1], meta=TRUE)
testthat::expect_equal(
  infer_format_ifti(readNifti(rs_nii_fnames[1])),
  c("nifti", "niftiImage")
)
testthat::expect_equal(
  infer_format_ifti(readNIfTI(rs_nii_fnames[1])),
  c("nifti", "nifti")
)
testthat::expect_equal(
  infer_format_ifti(array(rnorm(3*4*5*6), dim=c(3,4,5,6))),
  c("nifti", NA)
)

# data
testthat::expect_equal(
  infer_format_ifti(matrix(rnorm(15), nrow=5, ncol=3)),
  c("data", NA)
)
z <- as.matrix_ifti(matrix(rnorm(15), nrow=5, ncol=3), meta=TRUE)

# bad cases
testthat::expect_warning(infer_format_ifti(rs_cii_fnames))
testthat::expect_warning(infer_format_ifti("notAGoodFile.ii"))
testthat::expect_warning(infer_format_ifti(5))
