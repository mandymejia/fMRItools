#' Expand realignment parameters (RPs)
#' 
#' Compute the squares, differences, and square differences of each RP 
#'  timeseries.
#' 
#' @param RPs A \eqn{T \times N} numeric matrix, where \eqn{T} is the number of
#'  timepoints and \eqn{N} is the number of RPs (typically six) to expand.
#' @return A \eqn{T \times 4N} numeric matrix, with the first \eqn{N} columns
#'  being the original \code{RPs}, the next \eqn{N} being the differences,
#'  the next \eqn{N} being the squares, and the last \eqn{N} being the
#'  squared differences.
#' 
#' @export
#' 
expand_RPs <- function(RPs) {
  stopifnot(is.matrix(RPs))
  stopifnot(is.numeric(RPs))
  diffs <- rbind(0, diff(RPs))
  cbind(RPs, diffs, RPs^2, diffs^2)
}
