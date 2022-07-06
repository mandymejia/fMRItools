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
#' @return The input matrix with its columns centered and scaled.
#'
#' @importFrom robustbase rowMedians
scale_med <- function(mat, TOL=1e-8, drop_const=TRUE){
  # Transpose.
  mat <- t(mat)

  #	Center.
  mat <- mat - c(rowMedians(mat, na.rm=TRUE))

  # Scale.
  mad <- 1.4826 * rowMedians(abs(mat), na.rm=TRUE)
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