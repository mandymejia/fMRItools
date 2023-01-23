#' Format a path
#'
#' Normalize and validate a path (optionally, within a certain directory).
#'
#' @param path The path to normalize.
#' @param dir (Optional) the directory to append to the beginning of the path.
#'  \code{NULL} (default) to not append any directory, leaving \code{path}
#'  unchanged.
#' @param mode The mode for \code{\link{file.access}} to verify existence,
#'  writing permission, or reading permission. Use NA (default) to not perform
#'  any is.
#'
#' @return The normalized path, or \code{NULL} if the path was \code{NULL}.
#'
#' @keywords internal
#' 
format_path <- function(path, dir=NULL, mode=NA) {

  # Do nothing if the path is NULL.
  if (is.null(path)) { return(path) }

  # Append dir if provided.
  if (!is.null(dir)) { path <- file.path(dir, path) }
  path <- normalizePath(path, mustWork=FALSE)

  # Get the full file path (for Linux: previous normalizePath() does not get
  #   full file path if dir did not exist.)
  path <- file.path(
    normalizePath(dirname(path), mustWork=FALSE),
    basename(path)
  )

  # Check existence/writing permission/reading permission of the path.
  #   [NOTE]: This goes against this advice: 
  #   "Please note that it is not a good idea to use this
  #   function to test before trying to open a file. On a multi-tasking system,
  #   it is possible that the accessibility of a file will change between the
  #   time you call file.access() and the time you try to open the file. It is
  #   better to wrap file open attempts in try.
  stopifnot(all(mode %in% c(NA, 0, 2, 4)))
  for(m in mode) {
    if (is.na(mode)) { next }
    if (any(file.access(dirname(path), m) != 0)) {
      stop(paste0(
        "The directory \"", dirname(path), "\"",
        c(
          " doesn't exist. ", "",
          " is not writeable. Does it exist? ", "",
          "is not readable. Does it exist? "
        )[m+1],
        "Check and try again.\n"
      ))
    }
  }

  path
}

#' Is this an existing file path?
#'
#' Simple check if something is an existing file.
#'
#' @param x The potential file name
#'
#' @return Logical. Is \code{x} an existing file?
#'
#' @keywords internal
#'
is.fname <- function(x){
  is_1(x, "character") && file.exists(x) && !dir.exists(x)
}

#' Format a path for \code{\link{system}}
#' 
#' Right now, it uses \code{shQuote}
#'
#' @param R_path The name of the file. It should be properly formatted: if it
#'  exists, \code{file.exists(R_path)} should be \code{TRUE}.
#'
#' @return The name of the file
#'
#' @keywords internal
#' 
sys_path <- function(R_path) {
  shQuote(path.expand(R_path))
}