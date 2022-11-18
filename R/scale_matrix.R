#' Robust scaling
#' 
#' Centers and scales the columns of a matrix robustly
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). Constant-valued columns are set to \code{NA}
#' (or removed if \code{drop_const}) and a warning is raised. If all 
#' MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#' @param TOL minimum MAD to consider a column non-constant.
#'  Default: \code{1e-8}
#' @param drop_const Drop
#'
#' @export
#' @return The input matrix with its columns centered and scaled.
scale_med <- function(mat, TOL=1e-8, drop_const=TRUE){
  # Transpose.
  mat <- t(mat)

  #	Center.
  mat <- mat - c(rowMedians2(mat, na.rm=TRUE))

  # Scale.
  mad <- 1.4826 * rowMedians2(abs(mat), na.rm=TRUE)
  mad <- as.numeric(mad)
  const_mask <- mad < TOL
  if (any(const_mask)) {
    if (all(const_mask)) {
    stop("All columns are zero-variance.\n")
    } else {
      warning(paste0(
        "Warning: ", sum(const_mask),
        " constant columns (out of ", length(const_mask),
        " ). These will be removed.\n"
      ))
    }
  }
  mad <- mad[!const_mask]
  mat[const_mask,] <- NA
  mat[!const_mask,] <- mat[!const_mask,] / mad

  # Revert transpose.
  mat <- t(mat)

  if (drop_const) { mat <- mat[!const_mask,] }

  mat
}

#' Scale the BOLD timeseries
#'
#' @param BOLD Input fMRI data (V x T)
#' @param scale Option for scaling units.
#'
#' 	If \code{"auto"} (default), will use mean scaling except if demeaned data
#' 	is detected, in which case sd scaling will be used instead.
#'
#' 	\code{"mean"} scaling will scale the data to percent local signal change.
#'
#' 	\code{"sd"} scaling will scale the data by local standard deviation.
#'
#' 	\code{"none"} will only center the data, not scale it.
#' @param transpose Check orientation of data, which, if TRUE, will transpose
#' 	the data when the number of time points is greater than the number of voxels.
#' 	Note: this is not always true for subcortical regions.
#'
#' @return Scale to units of percent local signal change and centers
#'
#' @importFrom stats var
#' @export
scale_timeseries <- function(BOLD, scale=c("auto", "mean", "sd", "none"), transpose = TRUE){

	BOLD <- as.matrix(BOLD)
	nvox <- nrow(BOLD)
	ntime <- ncol(BOLD)

	scale <- match.arg(scale, c("auto", "mean", "sd", "none"))

	# Check orientation, send warning message and transpose if necessary.
	if((ntime > nvox) & transpose){
		warning('More columns than rows. Transposing matrix so rows are data locations and columns are time points')
		BOLD <- t(BOLD)
		nvox <- nrow(BOLD)
		ntime <- ncol(BOLD)
	}

	# Get `v_means`, the mean over time for each location (the mean image)
	v_means <- rowMeans(BOLD, na.rm=TRUE)
	v_means_min <- min(v_means, na.rm = TRUE)

	# Determine `"auto"` scaling.
	if (scale == "auto") {
		scale <- if (v_means_min > 1) { "mean" } else { "sd" }
		# cat("Using", scale, "scaling.\n")
	}

	# Center and scale.
	BOLD <- BOLD - v_means
	if (scale == "mean") {
		if (v_means_min < .1) {
			stop("Some local means are less than 0.1. Please set `scale_BOLD = 'none'` or set `meanTol` > .1 in the call to BayesGLM.")
		} else if (v_means_min < 1) {
			warning("Scaling to percent signal change when locations have means less than 1 may cause errors or produce aberrant results.")
		}
		BOLD <- t(100*BOLD / v_means)
	} else if (scale == "sd") {
		v_sd <- sqrt(apply(BOLD, 1, var, na.rm=TRUE))
		v_sd[is.na(v_sd)] <- 0
		if (min(v_sd) < 1e-6) {
			stop("Some local sds are less than 1e-6. Please set `scale_BOLD = 'none' or set `varTol` > 1e-3 in the call to BayesGLM.")
		}
		BOLD <- t(BOLD / v_sd)
	} else {
		BOLD <- t(BOLD)
	}

	BOLD
}

#' Scale the design matrix
#'
#' @param design_mat The original (unscaled) design matrix that is T x K, where
#'     T is the number of time points, and k is the number of task covariates
#'
#' @return A scaled design matrix
#' 
#' @export
#' 
# @examples
# # Task 1
# t1 <-
#   specifydesign(
#     onsets = seq(0, 200, by = 40),
#     durations = 1,
#     totaltime = 200,
#     TR = 1,
#     effectsize = 1.3,
#     conv = "double-gamma",
#     param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
#   )
# # Task 2
# t2 <-
#   specifydesign(
#     onsets = seq(20, 200, by = 40),
#     durations = 1,
#     totaltime = 200,
#     TR = 1,
#     effectsize = 1.3,
#     conv = "double-gamma",
#     param = list(list(a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, c = 0.15))
#   )
# A <- cbind(t1,t2) # This is the design matrix
# B <- scale_design_mat(A)
scale_design_mat <- function(design_mat) {
  if (!inherits(design_mat, "matrix")) stop("The design matrix must be a matrix class object.")
  output_mat <- apply(design_mat,2,function(task) {
    returned_col <- task / max(task)
    returned_col <- returned_col - mean(returned_col)
    return(returned_col)
  })
  return(output_mat)
}