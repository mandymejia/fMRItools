#' Match user inputs to expected values
#'
#' Match each user input to an expected/allowed value. Raise a warning if either
#'  several user inputs match the same expected value, or at least one could not
#'  be matched to any expected value. \code{ciftiTools} uses this function to
#'  match keyword arguments for a function call. Another use is to match
#'  brainstructure labels ("left", "right", or "subcortical").
#'
#' @param user Character vector of user input. These will be matched to
#'  \code{expected} using \code{\link{match.arg}}.
#' @param expected Character vector of expected/allowed values.
#' @param fail_action If any value in \code{user} could not be
#'  matched, or repeated matches occurred, what should happen? Possible values
#'  are \code{"stop"} (default; raises an error), \code{"warning"}, and
#'  \code{"nothing"}.
#' @param user_value_label How to refer to the user input in a stop or warning
#'  message. If \code{NULL}, no label is used.
#'
#' @return The matched user inputs.
#'
#' @keywords internal
#'
match_input <- function(
  user, expected,
  fail_action=c("stop", "warning", "message", "nothing"),
  user_value_label=NULL) {

  fail_action <- match.arg(
    fail_action,
    c("stop", "warning", "message", "nothing")
  )
  unrecognized_FUN <- switch(fail_action,
                             stop=stop,
                             warning=warning,
                             message=message,
                             nothing=invisible
  )

  if (!is.null(user_value_label)) {
    user_value_label <- paste0("\"", user_value_label, "\" ")
  }
  msg <- paste0(
    "The user-input values ", user_value_label,
    "did not match their expected values. ",
    "Either several matched the same value, ",
    "or at least one did not match any.\n\n",
    "The user inputs were:\n",
    "\t\"", paste0(user, collapse="\", \""), "\".\n",
    "The expected values were:\n",
    "\t\"", paste0(expected, collapse="\", \""), "\".\n"
  )

  tryCatch(
    {
      matched <- match.arg(user, expected, several.ok=TRUE)
      if (length(matched) != length(user)) { stop() }
      return(matched)
    },
    error = function(e) {
      unrecognized_FUN(msg)
    },
    finally = {
      NULL
    }
  )

  invisible(NULL)
}

#' Do these character vectors match exactly?
#' 
#' Checks if a user-defined character vector matches an expected character
#'  vector. That is, they share the same lengths and entries in the same order.
#'  For vectors of the same lengths, the result is \code{all(a == b)}.
#' 
#' Attributes are ignored.
#'
#' @param user Character vector of user input. 
#' @param expected Character vector of expected/allowed values.
#' @param fail_action If any value in \code{user} could not be
#'  matched, or repeated matches occurred, what should happen? Possible values
#'  are \code{"message"} (default), \code{"warning"}, \code{"stop"}, and
#'  \code{"nothing"}.
#'
#' @return Logical. Do \code{user} and \code{expected} match?
#' 
#' @keywords internal
#' 
match_exactly <- function(
  user, expected,
  fail_action=c("message", "warning", "stop", "nothing")) {

  fail_action <- match.arg(fail_action, c("message", "warning", "stop", "nothing"))
  unrecognized_FUN <- switch(fail_action,
    message=message,
    warning=warning,
    stop=stop,
    nothing=invisible
  )

  msg <- paste0(
    "Mismatch between:\n",
    "\t\"", paste0(user, collapse="\", \""), "\".\n",
    "and:\n",
    "\t\"", paste0(expected, collapse="\", \""), "\".\n"
  )

  tryCatch(
    {
      if (length(user) != length(expected)) { stop("Different lengths.") }
      if (!all(user == expected)) { stop("At least one different entry.") }
    },
    error = function(e) {
      unrecognized_FUN(msg)
      return(FALSE)
    },
    finally = {}
  )

  return(TRUE)
}

#' Is this object the expected data type, and length one?
#' 
#' @param x The value to check
#' @param dtype The data type. Default: \code{"numeric"}. Also can be 
#'  \code{"logical"} or \code{"character"}
#' @return \code{TRUE} if \code{x} is \code{dtype} and length one. 
#' @keywords internal
#'
is_1 <- function(x, dtype=c("numeric", "logical", "character")){
  dtype <- match.arg(dtype, c("numeric", "logical", "character"))
  xFUN <- switch(dtype, 
    numeric=is.numeric, 
    logical=is.logical, 
    character=is.character
  )
  xFUN(x) && (length(x)==1)
}

#' Is this object a positive number? (Or non-negative)
#' 
#' @param x The value to check
#' @param zero_ok Is a value of zero ok?
#' @return Logical indicating if \code{x} is a single positive or non-negative 
#'  number
#' @keywords internal
#' 
is_posNum <- function(x, zero_ok=FALSE){
  is_1(x) && ((x>0) || (x==0 && zero_ok))
}

#' All integers?
#'
#' Check if a data vector or matrix is all integers.
#'
#' @param x The data vector or matrix
#' @keywords internal
#'
#' @return Logical. Is \code{x} all integers?
#'
all_integers <- function(x){
  if (!is.numeric(x)) { return(FALSE) }
  non_integer <- max(abs(x - round(x)))
  non_integer==0 && !is.na(non_integer)
}

#' Is this an integer?
#' 
#' @param x The putative integer
#' @param nneg Require \code{x>=0} (non-negative) too?
#' @return Logical indicating whether \code{x} is an integer
#' 
#' @keywords internal
is_integer <- function(x, nneg=FALSE){
  out <- FALSE
  if (is_1(x)) {
    if (x%%1==0) {
      if (x>=0 || !nneg) { out <- TRUE }
    }
  } 
  out
}

#' Is this numeric vector constant?
#' 
#' @param x The numeric vector
#' @param TOL minimum range of \code{x} to be considered non-constant.
#'  Default: \code{1e-8}
#' 
#' @return Is \code{x} constant? 
#' 
#' @keywords internal
is_constant <- function(x, TOL=1e-8) {
  stopifnot(is.numeric(x))
  abs(max(x) - min(x)) < TOL
}
