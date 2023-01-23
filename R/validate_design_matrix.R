#' Validate design matrix
#' 
#' Coerces \code{design} to a numeric matrix, and optionally checks that the
#'  number of rows is as expected. Sets constant-valued columns to 1, and scales
#'  all other columns. 
#' 
#' @param design The design matrix
#' @param T_ the expected number of rows in \code{design}. Default:
#'  \code{NULL} (no expected value to validate).
#' 
#' @return The (modified) design matrix
#' 
#' @export
validate_design_matrix <- function(design, T_=NULL) {
  class(design) <- "numeric"
  if (identical(design, 1)) { design <- matrix(1, nrow=T_) }
  design <- as.matrix(design)
  if (!is.null(T_)) { stopifnot(nrow(design) == T_) }
  # Set constant columns (intercept regressor) to 1, and scale the other columns.
  design_const_mask <- apply(design, 2, is_constant)
  design[,design_const_mask] <- 1
  design[,!design_const_mask] <- scale(design[,!design_const_mask])
  design
}