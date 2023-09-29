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
my_data <- "/Users/ddpham/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/O365-BL-STAT-StatMIND-Projects - General/Data"

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

gica_gii_fnames <- c(
  "/Users/ddpham/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/O365-BL-STAT-StatMIND-Projects - General/Data/MSC-OneSubject/melodic_IC.sep.L.func.gii",
  "/Users/ddpham/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/O365-BL-STAT-StatMIND-Projects - General/Data/MSC-OneSubject/melodic_IC.sep.R.func.gii"
)

# Do ---------------------------------------------------------------------------
q <- harmonize(rs_cii_fnames[seq(3)], gica_cii_fname)

z <- read_cifti(rs_cii_fnames[1], brainstructures="all")
z <- newdata_xifti(z, t(q$DR[[1]]$S))
plot(z)

# Do ---------------------------------------------------------------------------
q <- harmonize(
  list(list(rs_giiL_fnames[1], rs_giiR_fnames[1]), list(rs_giiL_fnames[2], rs_giiR_fnames[2])),
  gica_gii_fnames
)
