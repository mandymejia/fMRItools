#' Check design matrix
#' 
#' @param design The design matrix
#' 
#' @return The (modified) design matrix
#' 
#' @keywords internal
check_design_matrix <- function(design, T_=nrow(design)) {
  class(design) <- "numeric"
  if (identical(design, 1)) { design <- matrix(1, nrow=T_) }
  design <- as.matrix(design)
  stopifnot(nrow(design) == T_)
  # Set constant columns (intercept regressor) to 1, and scale the other columns.
  design_const_mask <- apply(design, 2, is_constant)
  design[,design_const_mask] <- 1
  design[,!design_const_mask] <- scale(design[,!design_const_mask])
  design
}