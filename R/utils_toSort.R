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

#' Create a mask based on vertices that are invalid
#'
#' @param BOLD A \eqn{V \times T} numeric matrix. Each row is a location.
#' @param meanTol,varTol Tolerance for mean and variance of each data location. 
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{-Inf} for \code{meanTol} (ignore), and \code{1e-6} for 
#'  {varTol}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats rowVars
#' @return A logical vector indicating valid vertices
#'
#' @keywords internal
make_mask <- function(BOLD, meanTol=-Inf, varTol=1e-6, verbose=TRUE){
  stopifnot(is.matrix(BOLD))

  mask_na <- mask_mean <- mask_var <- rep(TRUE, nrow(BOLD))
  # Mark columns with any NA or NaN values for removal.
  mask_na[apply(is.na(BOLD) | is.nan(BOLD), 1, any)] <- FALSE
  # Calculate means and variances of columns, except those with any NA or NaN.
  # Mark columns with mean/var falling under the thresholds for removal.
  mask_mean[mask_na][rowMeans(BOLD[mask_na,,drop=FALSE]) < meanTol] <- FALSE
  if (varTol > 0) {
    mask_var[mask_na][matrixStats::rowVars(BOLD[mask_na,,drop=FALSE]) < varTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA/NaN values.\n"))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean.\n"))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance.\n"))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}

#' Positive skew?
#'
#' Does the vector have a positive skew?
#'
#' @param x The numeric vector for which to calculate the skew. Can also be a matrix,
#'  in which case the skew of each column will be calculated.
#' @return \code{TRUE} if the skew is positive or zero. \code{FALSE} if the skew is negative.
#' @keywords internal
#'
#' @importFrom stats median
skew_pos <- function(x){
  x <- as.matrix(x)
  apply(x, 2, median, na.rm=TRUE) <= colMeans(x, na.rm=TRUE)
}

#' Sign match ICA results
#'
#' Flips all source signal estimates (S) to positive skew
#'
#' @param x The ICA results with entries \code{S} and \code{M}
#' @return \code{x} but with positive skew source signals
#' @keywords internal
#'
sign_flip <- function(x){
  stopifnot(is.list(x))
  stopifnot(("S" %in% names(x)) & ("M" %in% names(x)))
  spos <- skew_pos(x$M)
  x$M[,!spos] <- -x$M[,!spos]
  x$S[,!spos] <- x$S[,!spos]
  x
}

#' Center cols
#'
#' Efficiently center columns of a matrix. (Faster than \code{scale})
#'
#' @param X The data matrix. Its columns will be centered
#' @return The centered data
#' @keywords internal
colCenter <- function(X) {
  X - rep(colMeans(X), rep.int(nrow(X), ncol(X)))
}

#' Check \code{Q2_max}
#'
#' Check \code{Q2_max} and set it if \code{NULL}.
#'
#' @param Q2_max,nQ,nT The args
#' @return \code{Q2_max}, clamped to acceptable range of values.
#' @keywords internal
Q2_max_check <- function(Q2_max, nQ, nT){
  if (!is.null(Q2_max)) {
    if (round(Q2_max) != Q2_max || Q2_max <= 0) {
      stop('`Q2_max` must be `NULL` or a non-negative integer.')
    }
  } else {
    Q2_max <- pmax(round(nT*.50 - nQ), 1)
  }

  # This is to avoid the area of the pesel objective function that spikes close
  #   to rank(X), which often leads to nPC close to rank(X)
  if (Q2_max > round(nT*.75 - nQ)) {
    warning('`Q2_max` too high, setting to 75% of T.')
    Q2_max <- round(nT*.75 - nQ)
  }

  Q2_max
}

#' Unmask a matrix
#'
#' @param dat The data
#' @param mask The mask
#' @keywords internal
unmask_mat <- function(dat, mask){
  stopifnot(nrow(dat) == sum(mask))
  mdat <- matrix(NA, nrow=length(mask), ncol=ncol(dat))
  mdat[mask,] <- dat
  mdat
}

#' Undo a volumetric mask
#' 
#' Un-applies a mask to vectorized data to yield its volumetric representation.
#'  The mask and data should have compatible dimensions: the number of rows in
#'  \code{dat} should equal the number of locations within the \code{mask}.
#'  This is used for the subcortical CIFTI data.
#' 
#' @param dat Data matrix with locations along the rows and measurements along 
#'  the columns. If only one set of measurements were made, this may be a 
#'  vector.
#' @param mask Volumetric binary mask. \code{TRUE} indicates voxels inside the
#'  mask.
#' @param fill The value for locations outside the mask. Default: \code{NA}.
#'
#' @return The 3D or 4D unflattened volume array
#'
#' @export
#' 
unmask_subcortex <- function(dat, mask, fill=NA) {

  # Check that dat is a vector or matrix.
  if (is.vector(dat) || is.factor(dat)) { dat <- matrix(dat, ncol=1) }
  stopifnot(length(dim(dat)) == 2)

  # Check that mask is numeric {0, 1} or logical, and is 3D.
  if (is.numeric(mask)) {
    mask_vals <- unique(as.vector(mask))
    stopifnot(length(mask_vals) <= 2)
    stopifnot(all(mask_vals %in% c(0,1)))
    mask <- array(as.logical(mask), dim=dim(mask))
  }
  stopifnot(length(dim(mask)) == 3)

  # Other checks.
  stopifnot(is.vector(fill) && length(fill)==1)
  stopifnot(sum(mask) == nrow(dat))

  # Make volume and fill.
  vol <- array(fill, dim=c(dim(mask), ncol(dat)))
  for (ii in seq_len(ncol(dat))) {
    vol[,,,ii][mask] <- dat[,ii]
  }
  if (ncol(dat)==1) { vol <- vol[,,,1] }

  vol
}

#' Pad a 3D Array
#'
#' Pad a 3D array by a certain amount in each direction, along each dimension.
#'  This effectively undoes a crop.
#'
#' @param x A 3D array, e.g. \code{unmask_subcortex(xifti$data$subcort, xifti$meta$subcort$mask)}.
#' @param padding A \eqn{d \times 2} matrix indicating the number of 
#'  slices to add at the beginning and end of each of the d dimensions, e.g.
#'  \code{xifti$meta$subcort$mask_padding}.
#' @param fill Values to pad with. Default: \code{NA}.
#'
#' @return The padded array
#'
#' @keywords internal
#' 
pad_vol <- function(x, padding, fill=NA){
  stopifnot(length(dim(x))==3)
  new_dim <- vector("numeric", 3)
  for (ii in seq_len(3)) {
    new_dim[ii] <- dim(x)[ii] + padding[ii,1] + padding[ii,2]
  }
  y <- array(fill, dim=new_dim)
  y[
    seq(padding[1,1]+1, padding[1,1]+dim(x)[1]),
    seq(padding[2,1]+1, padding[2,1]+dim(x)[2]),
    seq(padding[3,1]+1, padding[3,1]+dim(x)[3])
  ] <- x
  y
}

#' @rdname pad_vol
#' 
uncrop_vol <- function(x, padding, fill=NA){
  pad_vol(x, padding, fill)
}

#' Convert coordinate list to volume
#' 
#' Converts a sparse coordinate list to its non-sparse volumetric representation.
#' 
#' @param coords The sparse coordinate list. Should be a data.frame or matrix
#'  with voxels along the rows and three or four columns. The first three 
#'  columns should be integers indicating the spatial coordinates of the voxel.
#'  If the fourth column is present, it will be the value used for that voxel.
#'  If it is absent, the value will be \code{TRUE} or \code{1} if \code{fill}
#'  is not those values, and \code{FALSE} or \code{0} if \code{fill} is. The
#'  data type will be the same as that of \code{fill}.
#'  \code{fill}. The fourth column must be logical or numeric.
#' @param fill Fill value for the volume. Must be logical or numeric. Default: 
#'  \code{FALSE}.
#' 
#' @return The volumetric data
#'
#' @keywords internal
#' 
coordlist_to_vol <- function(coords, fill=FALSE){
  stopifnot(length(fill)==1)
  if (is.logical(fill)) {
    logical_vals <- TRUE
  } else if (is.numeric(fill)) {
    logical_vals <- FALSE
  } else { stop("Fill value must be logical or numeric.") }

  stopifnot(is.matrix(coords) || is.data.frame(coords))
  stopifnot(dim(coords)[2] %in% c(3,4))
  if (dim(coords)[2] == 3) {
    val <- ifelse(
      as.numeric(fill) != 1,
      ifelse(logical_vals, TRUE, 1),
      ifelse(logical_vals, FALSE, 0)
    )
    coords <- cbind(coords, val)
  } else {
    if (any(coords[,4] == fill, na.rm=TRUE)) { 
      warning("The fill value occurs in the data.") 
    }
  }

  vol <- array(fill, dim=apply(coords[,1:3], 2, max, na.rm=TRUE))
  vol[as.matrix(coords[,1:3])] <- coords[,4]
  vol
}

#' Crop a 3D array
#' 
#' Remove empty (zero-valued) edge slices from a 3D array.
#'
#' @param x The 3D array to crop.
#'
#' @keywords internal
#' 
crop_vol <- function(x) {
  d <- length(dim(x))

  if (all(unique(as.vector(x)) %in% c(NA, 0))) { stop("The array is empty.") }

  padding <- matrix(NA, nrow=d, ncol=2)
  rownames(padding) <- strsplit(rep("ijklmnopqrstuvwxyz", ceiling(d/15)), "")[[1]][1:d]
  empty_slice <- vector("list", d)
  for (ii in 1:d) {
    empty_slice[[ii]] <- apply(x, ii, sum, na.rm=TRUE) == 0
    first_slice <- min(which(!empty_slice[[ii]]))
    last_slice <- max(which(!empty_slice[[ii]]))
    padding[ii,1] <- ifelse(
      first_slice != 1, 
      first_slice - 1, 
      0
    )
    padding[ii,2] <- ifelse(
      last_slice != length(empty_slice[[ii]]), 
      length(empty_slice[[ii]]) - last_slice, 
      0
    )
  }
  x <- x[!empty_slice[[1]], !empty_slice[[2]], !empty_slice[[3]]]

  return(list(data=x, padding=padding))
}

#' Get spatial locations of each voxel
#' 
#' Use subcortical metadata (mask, transformation matrix and units) to get
#'  voxel locations in 3D space.
#' 
#' @param mask,trans_mat,trans_units The subcortical metadata
#' @return A list: \code{coords} and \code{units}
#' 
#' @keywords internal
#' 
vox_locations <- function(mask, trans_mat, trans_units=NULL){
  list(
    coords = (cbind(which(mask, arr.ind=TRUE), 1) %*% trans_mat)[,seq(3)],
    trans_units = NULL
  )
}