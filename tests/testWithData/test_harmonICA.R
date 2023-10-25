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
dir_StatMIND <- file.path("/Users/ddpham/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/O365-BL-STAT-StatMIND-Projects - General")
my_data <- file.path(dir_StatMIND, "Data")

subjects <- c(100307, 100408, 100610)
rs_cii_fnames <- file.path(my_data, "templateICAr-smallDataset", c(
  paste0(subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
))
rs_nii_fnames <- file.path(my_data, "templateICAr-smallDataset", c(
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

gica_cii_fname <- file.path(
  my_data, "templateICAr-smallDataset",
  "melodic_IC.4k.dscalar.nii"
)

gica_gii_fnames <- file.path(dir_StatMIND, c(
  "Data/MSC-OneSubject/melodic_IC.sep.L.func.gii",
  "Data/MSC-OneSubject/melodic_IC.sep.R.func.gii"
))

# CIFTI ------------------------------------------------------------------------
q <- harmonize(rs_cii_fnames[seq(3)], brainstructures=c("left", "right"), gica_cii_fname, TR=.72, do_harmonize=FALSE)
y <- read_cifti(rs_cii_fnames[1], brainstructures=c("left", "right"))
z <- newdata_xifti(y, t(q$S[1,,]))
plot(y); plot(z)
rgl::close3d(); rgl::close3d()

# GIFTI2 -----------------------------------------------------------------------
rs_cii_sep <- lapply(as.list(rs_cii_fnames), separate_cifti, write_dir=tempdir())
rs_cii_sep_dat <- lapply(rs_cii_sep, function(q){as.list(q[c("cortexL", "cortexR")])})
rs_cii_sep_mwall <- as.list(rs_cii_sep[[1]][c("ROIcortexL", "ROIcortexR")])
gica_cii_sep <- separate_cifti(gica_cii_fname, write_dir=tempdir())
q2 <- harmonize(
  rs_cii_sep_dat[seq(3)],
  as.list(gica_cii_sep[c("cortexL","cortexR")]),
  rs_cii_sep_mwall, TR=.72, do_harmonize=FALSE
)
testthat::expect_equal(q, q2)

# GIFTI ------------------------------------------------------------------------
q3 <- harmonize(
  vapply(rs_cii_sep_dat[seq(3)], function(q){q[[1]]}, ''),
  gica_cii_sep["cortexL"],
  rs_cii_sep_mwall[["ROIcortexL"]],
  inds=seq(4), scale="local", scale_sm_FWHM=2, hpf=0, GSR=TRUE,
  TR=.72, do_harmonize=FALSE
)

cor(q3$S[1,3,], q2$S[1,3,seq(3670)])
