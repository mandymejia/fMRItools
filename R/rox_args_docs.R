#' fMRI data for \code{scrub} and \code{CompCor}
#' 
#' @param X Wide numeric data matrix (\eqn{T observations} by \eqn{V variables}, \eqn{T << V}).
#'  For example, if \code{X} represents an fMRI run, \eqn{T} should be the number
#'  of timepoints and \eqn{V} should be the number of brainordinate vertices/voxels.
#' 
#'  Or, a 4D array or NIFTI or file path to a NIFTI (\eqn{I} by \eqn{J} by \eqn{K} by \eqn{T} 
#'  observations), in which case \code{ROI_data} must be provided. 
#'  (The vectorized data will be \eqn{T timepoints} by \eqn{V_{in-mask} voxels})
#' 
#'  Or, a \code{ciftiTools} \code{"xifti"} object or a file path to a CIFTI
#'  (The vectorized data will be \eqn{T timepoints} by \eqn{V_{left+right+sub} greyordinates}).
#' @param ROI_data Indicates the data ROI. Allowed arguments depend on \code{X}:
#' 
#'  If \code{X} is a matrix, this must be a length \eqn{V} logical vector, where
#'  the data ROI is indicated by \code{TRUE} values. If \code{"infer"} (default), all 
#'  columns of \code{X} will be included in the data ROI (\code{rep(TRUE, V)}).
#' 
#'  If \code{X} is an array or NIFTI, this must be either a vector of values
#'  to expect for out-of-mask voxels in \code{X}, or a (file path to a) 3D NIFTI.
#'  In the latter case, each of the volume dimensions should match the first
#'  three dimensions of \code{X}. Voxels in the data ROI should be indicated by
#'  \code{TRUE} and all other voxels by \code{FALSE}. If \code{"infer"} (default),
#'  will be set to \code{c(0, NA, NaN)} (include all voxels which are not constant
#'  \code{0}, \code{NA}, or \code{NaN}).
#' 
#'  If \code{X} is a \code{"xifti"} this must be the \code{brainstructures}
#'  argument to \code{ciftiTools::read_cifti}. If \code{"infer"} (default),
#'  \code{brainstructures} will be set to \code{"all"} (use both left and right
#'  cortex vertices, and subcortical voxels).
#'
#'  If \code{NULL}, the data ROI will be empty. This is useful for obtaining just
#'  the noise ROI, if the data and noise are located in separate files.
#' @param ROI_noise Indicates the noise ROIs for aCompCor. Should be a list where
#'  each entry corresponds to a distinct noise ROI. The names of the list should
#'  be the ROI names, e.g. \code{"white_matter"} and \code{"csf"}. The expected
#'  formats of the list entries depends on \code{X}:
#' 
#'  For all types of \code{X}, \code{ROI_noise} entries can be a matrix of noise
#'  ROI data. The matrix should have \eqn{T} rows, with each column being a
#'  data location's timeseries.
#'  
#'  If \code{X} is a matrix, entries can also indicate a noise ROI within \code{X}.
#'  These entries must be a length \eqn{V} logical vector with \code{TRUE} values 
#'  indicating locations in \code{X} within that noise ROI. Since the ROIs must 
#'  not overlap, the masks must be mutually exclusive with each other, and with 
#'  \code{ROI_data}. 
#' 
#'  If \code{X} is an array or NIFTI, entries can also indicate a noise ROI within \code{X}.
#'  These entries must be a logical array or (file path to) a 3D NIFTI with the
#'  same spatial dimensions as \code{X}, and with \code{TRUE} values indicating
#'  voxels inside the noise ROI. Since the ROIs must not overlap, the masks must
#'  be mutually exclusive with each other, and with \code{ROI_data}. 
#' 
#'  (If \code{X} is a \code{"xifti"}, entries must be data matrices, since no 
#'  greyordinate locations in \code{X} are appropriate noise ROIs).
#' @name data_CompCor_Params
#' @keywords internal
NULL
#' noise parameters for CompCor
#' @param noise_nPC The number of principal components to compute for each noise
#'  ROI. Alternatively, values between 0 and 1, in which case they will 
#'  represent the minimum proportion of variance explained by the PCs used for
#'  each noise ROI. The smallest number of PCs will be used to achieve this 
#'  proportion of variance explained. 
#' 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its
#'  elements will be recycled. Default: \code{5} (compute the top 5 PCs for 
#'  each noise ROI).
#' @param noise_erosion The number of voxel layers to erode the noise ROIs by. 
#'  Should be a list or numeric vector with the same length as \code{ROI_noise}. 
#'  It will be matched to each ROI based on the name of each entry, or if the 
#'  names are missing, the order of entries. If it is an unnamed vector, its 
#'  elements will be recycled. Default: \code{NULL}, which will use a value of
#'  0 (do not erode the noise ROIs). Note that noise erosion can only be
#'  performed if the noise ROIs are volumetric.
#' @name noise_Params
#' @keywords internal
NULL