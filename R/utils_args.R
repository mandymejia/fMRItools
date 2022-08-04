#' Get kwargs
#' 
#' Get the names of the arguments of a function as a character vector.
#'
#' @param fun The function to get the argument names for.
#'
#' @return The names of the arguments of \code{fun} as a character vector
#'
#' @keywords internal
#' 
get_kwargs <- function(fun) {
  kwargs <- names(as.list(args(fun)))
  kwargs <- kwargs[seq(length(kwargs)-1)] # last is empty
  kwargs
}

#' Merges two kwargs 
#' 
#' Merge two kwarg lists. If a kwarg is present in both lists but with different
#'  values, an error is raised.
#' @param kwargsA The first list of kwargs.
#' @param kwargsB The second list of kwargs. If duplicates are present, the default
#'  message recommends the user to remove the kwarg here in favor of placing the
#'  correct one in \code{kwargsA}.
#' @param labelA (Optional) Descriptor of \code{kwargsA} for error statement. Default "first kwarg(s)".
#' @param labelB (Optional) Descriptor of \code{kwargsB} for error statement. Default "second kwarg(s)".
#' @param extraMsg (Optional) Extra text for error statement. "\[DEFAULT\]" (default) will use this message:
#'  "Note that a kwarg only has to be provided to one of these. Place the correct value in the first
#'  location and remove the kwarg from the second location".
#'
#' @return A list with the union of \code{kwargsA} and \code{kwargsB}
#'
#' @keywords internal
#' 
merge_kwargs <- function(kwargsA, kwargsB,
  labelA="first kwarg(s)", labelB="second kwarg(s)",
  extraMsg="[DEFAULT]") {

  # Identify repeated kwargs.
  repeatedB_bool <- names(kwargsB) %in% names(kwargsA)
  repeated <- names(kwargsB)[repeatedB_bool]
  # Stop if any repeated kwargs differ.
  kwargs_mismatch <- !mapply(identical, kwargsA[repeated], kwargsB[repeated])
  if (sum(kwargs_mismatch) > 0) {
    if(identical(extraMsg, "[DEFAULT]")){
      extraMsg <- "Note that a kwarg only has to be provided to one of these. \
        Place the correct value in the first location and remove the kwarg \
        from the second location"
    }
    stop(paste0(
      "A keyword argument(s) was provided twice with different values. Here is the kwarg(s) in disagreement:\n",
      "The ", labelA, " was:\n",
      "\"", paste0(kwargsA[kwargs_mismatch], collapse="\", \""), "\".\n",
      "The ", labelB, " was:\n",
      "\"", paste0(kwargsB[kwargs_mismatch], collapse="\", \""), "\".\n",
      extraMsg
    ))
  }
  kwargs <- c(kwargsA, kwargsB[!repeatedB_bool])
}