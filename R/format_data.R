#' Count each voxel's neighbors with value \code{TRUE}
#'
#' For each voxel in a 3D logical array, count the number of immediate neighbors
#'  with value \code{TRUE}.
#'
#' @param arr The 3D logical array.
#' @param pad Pad value for edge.
#' 
#' @return An array with the same dimensions as \code{arr}. Each voxel value
#'  will be the number of its immediate neighbors (0 to 6) which are \code{TRUE}.
#' 
#' @keywords internal
neighbor_counts <- function(arr, pad=FALSE){
  stopifnot(length(dim(arr)) == 3)
  stopifnot(all(unique(arr) %in% c(TRUE, FALSE)))
  arrPad <- pad_vol(arr, matrix(1, 3, 2), fill=pad)
  dPad <- dim(arrPad)
  # Look up, down, left, right, forward, behind (not diagonally)
  arrPad[1:(dPad[1]-2),2:(dPad[2]-1),2:(dPad[3]-1)] +
    arrPad[3:(dPad[1]),2:(dPad[2]-1),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),1:(dPad[2]-2),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),3:(dPad[2]),2:(dPad[3]-1)] +
    arrPad[2:(dPad[1]-1),2:(dPad[2]-1),1:(dPad[3]-2)] +
    arrPad[2:(dPad[1]-1),2:(dPad[2]-1),3:(dPad[3])]
}

#' Erode volumetric mask
#'
#' Erode a volumetric mask by a certain number of voxel layers. For each layer,
#'  any in-mask voxel adjacent to at least one out-of-mask voxel is removed
#'  from the mask. 
#'
#' Diagonal voxels are not considered adjacent, i.e. the voxel at (0,0,0) is not
#'  adjacent to the voxel at (1,1,0) or (1,1,1), although it is adjacent to 
#'  (1,0,0).
#'
#' @param vol The volume to erode. Out-of-mask voxels should be indicated by a 
#'  value in \code{out_of_mask_val}.
#' @param n_erosion The number of layers to erode the mask by.
#' @param out_of_mask_val A voxel is not included in the mask if and only if its
#'  value is in this vector. The first value in this vector will be used to replace
#'  the eroded voxels. Default: \code{NA}.
#'
#' @return The eroded \code{vol}. It is the same as \code{vol} but with eroded
#'  voxels replaced with the value \code{out_of_mask_val[1]}.
#'
#' @export
erode_vol <- function(vol, n_erosion=1, out_of_mask_val=NA){
  stopifnot(is_integer(n_erosion, nneg=TRUE))
  stopifnot(length(dim(vol)) == 3)
  if (n_erosion==0) { return(vol) }
  mask <- !array(vol %in% out_of_mask_val, dim=dim(vol))
  for (ii in seq(n_erosion)) {
    to_erode <- mask & neighbor_counts(mask, pad=TRUE) < 6
    mask[to_erode] <- FALSE
    vol[to_erode] <- out_of_mask_val[1]
  }
  vol
}

#' Format data for pscrub and CompCor
#'
#' @param X description
#' @param ROI_data,ROI_noise description
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
#'  0 (do not erode the noise ROIs).
#'
#' @return A list with components "X", "X_noise", "ROI_data", and "ROI_noise"
#' 
#' @keywords internal
format_data <- function(X, ROI_data="infer", ROI_noise=NULL, noise_nPC=5, noise_erosion=NULL){

  # TO DO: explain how to use ROI_data

  # if X is a file, read it.
  if (is.character(X)) {
    if (endsWith(X, ".dtseries.nii") | endsWith(X, ".dscalar.nii")) {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read `X`. Please install it", call. = FALSE)
      }
      if (identical(ROI_data, "infer")) { ROI_data <- "all" }
      X <- ciftiTools::read_cifti(X, brainstructures=ROI_data)
    } else if (endsWith(X, ".nii") | endsWith(X, ".nii.gz")) {
      X <- read_nifti(X)
    } else {
      stop("The argument `X`, '", X, "', does not look like a CIFTI or NIFTI file name.")
    }
  } 

  # X must be a matrix, array, or "xifti"
  if (is.matrix(X)) {
    T_ <- nrow(X); V_ <- ncol(X)
    if (T_ > V_) {
      warning(
        "Data matrix has more rows than columns. Check that observations\
        are in rows and variables are in columns."
      )
    }
    X_type <- "vector"
  } else if (is.array(X)) {
    if (length(dim(X))==3) { X <- array(X, dim=c(dim(X)[1:2], 1, dim(X)[3])) }
    stopifnot(length(dim(X))==4)
    T_ <- dim(X)[4]
    X_type <- "volume"
  } else if (inherits(X, "xifti")) {
    xifti_meta <- X$meta
    X <- t(as.matrix(X))
    T_ <- nrow(X); V_ <- ncol(X)
    X_type <- "xifti"
  } else {
    stop("`X` must be a matrix, array, NIFTI, path to a NIFTI, CIFTI, or path to a CIFTI.")
  }

  if (T_ < 2) { stop("There are less than two timepoints.") }

  # ROI_data
  if (is.null(ROI_data)) {
    if (X_type == "volume") {
      ROI_data <- X[,,,1] > Inf
    } else {
      ROI_data <- rep(FALSE, V_)
    }
    ROI_data_infer <- FALSE

  } else {
    ROI_data_infer <- identical(ROI_data, "infer")
    if (X_type == "vector") {
      if (ROI_data_infer) { ROI_data <- rep(TRUE, V_) }
      ROI_data <- as.vector(ROI_data)
      if (length(ROI_data) != V_) { 
        stop("The `ROI_data` must be a logical vector with the same length as columns in the data.") 
      }
      ROI_data <- as.logical(ROI_data)
    } else if (X_type == "volume") {
      if (ROI_data_infer) { ROI_data <- c(0, NA, NaN) }
      if (is.character(ROI_data)) {
        ROI_data <- read_nifti(ROI_data)
      }
      if (is.vector(ROI_data)) {
        if (length(ROI_data) != length(unique(ROI_data))) {
          stop(
            "`X` is a volume, and `ROI_data` is a vector, so `ROI_data` should be\
            values which out-of-mask voxels take on. But, the values of `ROI_data`\
            are not unique."
          )
        }
        ROI_data <- apply(array(X %in% ROI_data, dim=dim(X)), 1:3, function(x){!all(x)})
      } else if (is.array(ROI_data)) {
        if (length(dim(ROI_data))==2) { ROI_data <- array(ROI_data, dim=c(dim(ROI_data), 1)) }
        if (all(dim(ROI_data) != dim(X)[1:3])) { 
          stop("The `ROI_data` must have the same dimensions as the first three dimensions of `X`.")
        }
        mode(ROI_data) <- "logical"
      }
    } else if (X_type == "xifti") {
      ROI_data <- rep(TRUE, V_)
    } else { stop("Internal error: unrecognized `X_type`") }
    if (sum(ROI_data) == 0) { warning("The data ROI was empty.\n") }
  }

  # ROI_noise
  if (!is.null(ROI_noise)) { #&& length(ROI_noise)>0) {
    if (!is.list(ROI_noise)) { ROI_noise <- list(Noise1=ROI_noise) }
    if (is.null(names(ROI_noise))) { names(ROI_noise) <- paste0("Noise", 1:length(ROI_noise)) }
    if (length(names(ROI_noise)) != length(unique(names(ROI_noise)))) {
      stop("The `ROI_noise` names must be unique.")
    }
    stopifnot(!any(names(ROI_noise) == "data"))

    # noise_nPC
    if (length(noise_nPC) == 1) {
      noise_nPC <- as.list(rep(noise_nPC, length(ROI_noise)))
    } else {
      noise_nPC <- as.list(noise_nPC)
    }
    names(noise_nPC) <- names(ROI_noise)

    if (is.null(names(noise_nPC))) {
      noise_nPC <- noise_nPC[rep(1:length(noise_nPC), length(ROI_noise))[1:length(ROI_noise)]]
      names(noise_nPC) <- names(ROI_noise)
    }
    else {
      if (all(sort(names(noise_nPC)) != sort(names(ROI_noise)))) {
        stop("The names of `noise_nPC` do not match those of `ROI_noise`.")
      }
    }  
    noise_nPC <- noise_nPC[names(ROI_noise)]

    # noise_erosion
    if (X_type == "volume") {
      if (is.null(noise_erosion)) { 
        noise_erosion = 0
      } else {
        if (!all(noise_erosion==0) & !any(sapply(ROI_noise, is.array))) {
          warning("`noise_erosion` was provided, but there are no array/NIFTI noise ROIs to erode.\n")
        }
      }
      noise_erosion <- as.list(noise_erosion)
      if (is.null(names(noise_erosion))) {
        noise_erosion <- noise_erosion[rep(1:length(noise_erosion), length(ROI_noise))[1:length(ROI_noise)]]
        names(noise_erosion) <- names(ROI_noise)
      } else {
        if (all(sort(names(noise_erosion)) != sort(names(ROI_noise)))) {
          stop("The names of `noise_erosion` do not match those of `ROI_noise`.")
        }
      }
      noise_erosion <- noise_erosion[names(ROI_noise)]
      
    } else {
      if (!is.null(noise_erosion)) { 
        warning(
          "Erosion requires volumetric data, but the data is not volumetric.\
          No erosion will happen.\n"
        ) 
      }
    }

    X_noise <- vector("list", length(ROI_noise)); names(X_noise) <- names(ROI_noise)
    for (ii in 1:length(ROI_noise)) {
      if (is.null(ROI_noise[[ii]])) { ROI_noise[ii] <- list(NULL); next }
      if (X_type == "vector") {
        if (is.vector(ROI_noise[[ii]])) {
          stopifnot(length(ROI_noise[[ii]]) == V_)
          ROI_noise[[ii]] <- as.logical(ROI_noise[[ii]])
          X_noise[[ii]] <- X[,ROI_noise[[ii]]]
        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[ii] <- list(NULL)
        } else {
          stop(
            "Each entry in `ROI_noise` must be a logical vector, or matrix\
            with the same number of rows as timepoints in `X`."
          )
        }
      } else if (X_type == "volume") {
        if (is.character(ROI_noise[[ii]])) {
          if (!file.exists(ROI_noise[[ii]])) { 
            stop(paste(
              "The `ROI_noise` entry", ROI_noise[[ii]], "is not an existing file."
            )) 
          }
          ROI_noise[[ii]] <- read_nifti(ROI_noise[[ii]])
        }
        if (is.matrix(ROI_noise[[ii]])) { 
          ROI_noise[[ii]] <- array(ROI_noise[[ii]], dim=c(dim(ROI_noise[[ii]]), 1))
        }
        if (is.array(ROI_noise[[ii]])) {
          stopifnot(all(dim(ROI_noise[[ii]]) == dim(X)[1:3]))
          ROI_noise[[ii]][,,] <- as.logical(ROI_noise[[ii]]) * 1
          ROI_noise[[ii]] <- erode_vol(ROI_noise[[ii]], noise_erosion[[ii]], c(-1, 0, NA, NaN))
          X_noise[[ii]] <- t(matrix(X[ROI_noise[[ii]] > 0], ncol=T_))

        } else if (is.matrix(ROI_noise[[ii]])) {
          stopifnot(nrow(ROI_noise[[ii]]) == T_)
          X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[ii] <- list(NULL)
        } else {
          stop(
            "Each entry in `ROI_noise` must be a logical array, or matrix\
            with the same number of rows as timepoints in `X`."
          ) 
        }
      } else if (X_type == "xifti") {
        stopifnot(is.matrix(ROI_noise[[ii]]))
        stopifnot(nrow(ROI_noise[[ii]]) == T_)
        X_noise[[ii]] <- ROI_noise[[ii]]; ROI_noise[ii] <- list(NULL)
      } else { stop("Internal error: unrecognized `X_type`") }
      if (ncol(X_noise[[ii]]) == 0) { 
        warning(paste("The noise ROI", names(ROI_noise)[ii], "is empty."))
      }
    }
    ROI_noise <- ROI_noise[!vapply(ROI_noise, is.null, FALSE)]

    if (!all(sapply(ROI_noise, is.null))) {
      # check that ROI are mutually exclusive
      all_ROI_noises <- apply(do.call(rbind, ROI_noise) > 0, 2, sum)
      if (!all(all_ROI_noises < 2)) {
        stop("The noise ROIs must all be mutually exclusive.")
      }
      all_ROI_noises <- all_ROI_noises > 0
      if (ROI_data_infer) { 
        ROI_data[all_ROI_noises] <- FALSE
      } else {
        if (any(all_ROI_noises & as.vector(ROI_data))) {
          warning("The noise ROIs overlapped with the data ROI. Labeling overlapped voxels as noise.")
          if (X_type == "volume") {
            ROI_data[,,] <- ROI_data & !all_ROI_noises
          } else {
            ROI_data <- ROI_data & !all_ROI_noises
          }
        }
      }
    }
    if (length(ROI_noise) == 0) { ROI_noise <- NULL }

  } else {
    X_noise <- NULL; noise_nPC <- NULL
  }

  # make sure nPCs greater than the rank of the noise ROIs!
  # similarly for data
  # ...

  if (X_type == "volume") {
    X <- t(matrix(X[ROI_data], ncol=dim(X)[4]))
  } else {
    X <- X[,ROI_data]
  }

  list(
    X=X, X_noise=X_noise, 
    ROI_data=ROI_data, ROI_noise=ROI_noise, 
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )
}

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
  if(!(length(x)==1 & is.character(x))){ return(FALSE) }
  file.exists(x) & !dir.exists(x)
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

#' Match user inputs to expected values
#'
#' Match each user input to an expected/allowed value. 
#' 
#' Raise a warning if either
#'  several user inputs match the same expected value, or at least one could not
#'  be matched to any expected value. \code{ciftiTools} uses this function to
#'  match keyword arguments for a function call. Another use is to match
#'  brainstructure labels ("left", "right", or "subcortical").
#'
#' @param user Character vector of user input. These will be matched to
#'  \code{expected} using \code{match.arg()}.
#' @param expected Character vector of expected/allowed values.
#' @param fail_action If any value in \code{user} could not be
#'  matched, or repeated matches occurred, what should happen? Possible values
#'  are \code{"stop"} (default; raises an error), \code{"warning"}, and
#'  \code{"nothing"}.
#' @param user_value_label How to refer to the user input in a stop or warning
#'  message. If \code{NULL}, no label is used.
#'
#' @return The matched user inputs
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

#' Identify Boundary Layers.
#'
#' Identify the vertices within \code{boundary_width} edges of the input mask.
#'  The mesh must be triangular.
#' 
#' @param faces a V x 3 matrix of integers. Each row defines a face by the index
#'  of three vertices.
#' @inheritParams mask_Param_vertices
#' @param boundary_width a positive integer representing the width of the boundary
#'  to compute. The furthest vertices from the input mask will be this number of
#'  edges away from the closest vertex in the input mask. Default: \code{10}.
#' 
#' @return a length-V numeric vector. Each entry corresponds to the vertex
#'  with the same index. For vertices within the boundary, the value will be the
#'  number of vertices away from the closest vertex in the input mask.
#'  Vertices inside the input mask but at the edge of it (touching vertices with
#'  value 1) will have value 0. All other vertices will have value -1.
#'
#' @keywords internal 
boundary_layers <- function(faces, mask, boundary_width=10){
  s <- ncol(faces)
  v <- max(faces)
  # For quads, boundary_layers() would count opposite vertices on a face as
  #   adjacent--that's probably not desired.
  stopifnot(s == 3)

  stopifnot(boundary_width > 0)

  verts_mask <- which(mask)
  b_layers <- rep(-1, v)
  for (ii in seq(1, boundary_width)) {
    # Identify vertices that share a face with in-mask vertices, or vertices
    #   in lower layers. Those faces occupy the next boundary layer.
    ## Vertices inside the mask, or a lower layer
    vert_maskorlower <- unique(c(verts_mask, which(b_layers > 0)))
    ## The number of "in-mask or lower-layer" vertices in each face
    face_n_maskorlower <- rowSums(matrix(faces %in% vert_maskorlower, ncol=s))
    ## Faces with a mix of "in-mask or lower-layer" vertices, and new vertices
    faces_adj_mask <- face_n_maskorlower > 0 & face_n_maskorlower < s
    if (!any(faces_adj_mask)) { break }
    ## Vertices in faces with a mix
    verts_adj <- unique(as.vector(faces[faces_adj_mask,]))

    ## For the first layer, in-mask vertices that are part of a mixed face
    ##  are "layer 0"
    if (ii == 1) {
      b_layers[verts_adj[verts_adj %in% which(mask)]] <- 0
    }

    ## The new vertices that are part of the mixed faces make up this layer.
    verts_adj <- verts_adj[!(verts_adj %in% vert_maskorlower)]
    b_layers[verts_adj] <- ii
  }

  b_layers
}

#' Vertex Adjacency Matrix
#'
#' Make adjacency matrix between two sets of vertices.
#'
#' @param faces a V x 3 matrix of integers. Each row defines a face by the index
#'  of three vertices.
#' @param v1,v2 The first and second set of vertices. These are logical vectors
#'  the same length as \code{vertices} indicating the vertices in each set.
#'  If \code{v2} is \code{NULL} (default), set \code{v2} to \code{v1}. Can
#'  alternatively be a vector if integers corresponding to vertex indices.
#'
#' @return Adjacency matrix
#' 
#' @keywords internal 
vert_adjacency <- function(faces, v1, v2=NULL){
  v_all <- unique(as.vector(faces))
  # Arguments.
  if (is.logical(v1)) { v1 <- which(v1) }
  v1 <- sort(unique(v1))
  stopifnot(all(v1 %in% v_all))
  if (is.null(v2)) {
    v2 <- v1
  } else {
    if (is.logical(v2)) { v2 <- which(v2) }
    v2 <- sort(unique(v2))
    stopifnot(all(v2 %in% v_all))
  }

  # # Faces with at least one vertex in v1, and at least one vertex in v2.
  # #   This step reduces the number of faces we check in the next step.
  # #   It also used the subset index.
  # f_btwn <- apply(matrix(faces %in% v1, ncol=3), 1, any)
  # f_btwn_mask <- f_btwn & apply(matrix(faces %in% v2, ncol=3), 1, any)
  # f_btwn <- faces[f_btwn_mask,]
  v_reidx <- rep(NA, max(as.numeric(faces)))
  v_reidx[v1] <- 1:length(v1)
  v_reidx[v2] <- 1:length(v2)
  f_reidx <- matrix(v_reidx[as.vector(faces)], ncol=ncol(faces))

  # Check each pair of vertices in each face.
  adj <- matrix(FALSE, nrow=length(v1), ncol=length(v2))
  for (ii in 1:3) {
    for (jj in 1:3) {
      if (ii == jj) { next }
      # Mark adjacency betwen (v1, v2) pairs sharing a face.
      v_pairs <- f_reidx[(faces[,ii] %in% v1) & (faces[,jj] %in% v2), c(ii,jj)]
      adj[v_pairs] <- TRUE
    }
  }

  # Add "v" to row/colnames to not confuse with numeric index.
  rownames(adj) <- paste("v", v1); colnames(adj) <- paste("v", v2)
  adj
}

#' Order Vertices on Circular Manifold
#'
#' Order vertices on circular manifold by radians (after 2D CMDS projection).
#'
#' @inheritParams vertices_Param
#' 
#' @return Index ordering of \code{vertices}
#' 
#' @importFrom stats cmdscale dist
#' @keywords internal 
radial_order <- function(vertices){
  # Use CMDS to project onto 2-dimensional subspace. Scale each dimension.
  x <- scale(cmdscale(dist(vertices)))
  # Remove zero-values to stay in the domain of the trig functions.
  x <- matrix(ifelse(abs(x) < 1e-8, sign(x)*1e-8, x), ncol=ncol(x))
  # Order by the radians counter-clockwise from positive x-axis.
  # https://stackoverflow.com/questions/37345185/r-converting-cartesian-to-polar-and-sorting
  order(order(ifelse(
    x[,1] < 0,
    atan(x[,2] / x[,1]) + pi,
    ifelse(x[,2] < 0 , atan(x[,2] / x[,1]) + 2*pi, atan(x[,2] / x[,1]))
  )))
}

#' Apply Mask With Boundary To Mesh
#'
#' Make a boundary around a mask with two levels of decimation, and apply to a mask.
#'
#' The boundary consists of a \code{width1}-vertex-wide middle region and a
#'  \code{width2}-vertex-wide outer region, for a total of \code{width1 + width2} layers
#'  of vertices surrounding the input mask. In the first layer, every \code{k1}
#'  vertex within every \code{k1} layer (beginning with the innermost
#'  layer) is retained; the rest are discarded. In the second layer, every
#'  \code{k2} vertex within every \code{k2} layer (beginning with the innermost
#'  layer) is retained; the rest are discarded. It is recommended to make \code{width1}
#'  a multiple of \code{k1} and \code{width2} a multiple of \code{k2}.
#'
#' Default boundary: a 4-vertex wide middle region with triangles twice as long,
#'  and a 6-vertex wide outer region with triangles three times as long.
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mask_Param_vertices
#' @param width1,width2 the width of the middle/outer region. All vertices in the middle/outer region
#'  are between 1 and \code{width1} edges away from the closest vertex in \code{mask}/middle region.
#' @param k1,k2 roughly, the triangle size multiplier. Every \code{k1}/\code{k2} vertex within
#'  every \code{k1}/\code{k2} layer (beginning with the innermost layer) will be retained;
#'  the rest will be discarded. If the mesh originally has triangles of regular
#'  size, the sides of the triangles in the middle/outer region will be about
#'  \code{k1}/\code{k2} as long.
#'
#' @return A new mesh (list with components vertices and faces)
#' 
#' @importFrom stats quantile
#'
#' @keywords internal 
mask_with_boundary <- function(vertices, faces, mask, width1=4, k1=2, width2=6, k2=3){

  # ----------------------------------------------------------------------------
  # Check arguments ------------------------------------------------------------
  # ----------------------------------------------------------------------------
  width <- width1 + width2
  V_ <- nrow(vertices)
  F_ <- nrow(faces)
  s <- ncol(faces)
  stopifnot(s==3)

  faces2 <- list(
    rm = rep(FALSE, F_),
    add = NULL
  )

  # ----------------------------------------------------------------------------
  # Pre-compute layers and vertex adjacency matrix between neighbor layers. ----
  # ----------------------------------------------------------------------------
  b_layers <- boundary_layers(faces, mask, width)
  b_adjies <- vector("list", width)
  for (ii in 1:length(b_adjies)){
    b_adjies[[ii]] <- vert_adjacency(
      faces,
      v1 = which(b_layers == ii-1),
      v2 = which(b_layers == ii)
    )
  }

  # ----------------------------------------------------------------------------
  # Working outward from the mask, collect info on each layer (and previous), --
  # ----------------------------------------------------------------------------
  lay_idxs <- c(0, seq(1, width1, k1))
  if (width2 != 0) { lay_idxs <- c(lay_idxs, seq(width1+1, width1+width2, k2)) }
  lay_k <- c(1, rep(k1, width1), rep(k2, width2))

  # At each iteration, we will calcuate these "layer facts":
  # Note: verts_pre must be in radial order: lay$verts[lay$rad_order],]
  get_layer_facts <- function(vertices, faces, b_layers, lay_idx, k, verts_pre=NULL){
    lay <- list(idx = lay_idx)
    ## The vertices in the layer
    lay$verts <- which(b_layers == lay$idx)
    ## The number of vertices
    lay$V1 <- length(lay$verts)
    ## Faces whose vertices are entirely in the layer
    lay$faces_complete <- apply(matrix(faces %in% lay$verts, ncol=s), 1, all)
    ## The radial ordering of the vertices in the layer
    lay$rad_order <- radial_order(vertices[lay$verts,])
    ## The vertices in radial order
    lay$verts_rad <- vertices[lay$verts[order(lay$rad_order)],]
    if (!is.null(verts_pre)) {
      ## Adjust radial ordering so first vertex in this layer is closest to
      ##  first vertex in pre layer.
      rad_first <- which.min(
        apply(t(lay$verts_rad) - verts_pre[1,], 2, norm)
      )
      lay$rad_order <- ((lay$rad_order - rad_first) %% length(lay$rad_order)) + 1
      lay$verts_rad <- lay$verts_rad[c(
        seq(rad_first, nrow(lay$verts_rad)),
        seq(1, rad_first-1)
      ),]
      ## Flip direction of radial ordering (clockwise/counter-clockwise, or rather
      ##  left/right to match the pre layer)
      even_sample <- function(X, s){ X[as.integer(floor(quantile(1:nrow(X), probs=seq(1,s)/s))),] }
      pre_samp <- even_sample(verts_pre, 12)
      lay_samp <- even_sample(lay$verts_rad, 12)
      flip <- mean(apply(pre_samp - lay_samp, 1, norm)) > mean(apply(pre_samp[nrow(pre_samp):1,] - lay_samp, 1, norm))
      if (flip) {
        lay$rad_order <- (length(lay$rad_order) - lay$rad_order + 2) %% length(lay$rad_order)
        lay$rad_order[lay$rad_order==0] <- length(lay$rad_order)
        lay$verts_rad <- lay$verts_rad[nrow(lay$verts_rad):1,]
      }
    }
    ## Only keep every kth
    lay$verts_rad <- lay$verts_rad[seq(1, nrow(lay$verts_rad), k),]
    lay$V2 <- nrow(lay$verts_rad)
    lay
  }

  lay_ii <- get_layer_facts(
    vertices, faces, b_layers,
    lay_idxs[1], lay_k[1],
  )

  for (ii in 1:(length(lay_idxs)-1)) {
    # Calculate layer facts for this layer (and keep the facts for the previous one).
    lay_pre <- lay_ii
    lay_ii <- get_layer_facts(
      vertices, faces, b_layers,
      lay_idxs[ii+1], lay_k[ii+1],
      lay_pre$verts_rad
    )
    # plot_3d(lay_ii$verts_rad, 1:nrow(lay_ii$verts_rad))
    # plot_3d(vertices[lay_ii$verts,], lay_ii$rad_order)

    # --------------------------------------------------------------------------
    # Remove faces between the layers. -----------------------------------------
    # --------------------------------------------------------------------------
    faces_btwn <- apply(
      matrix(faces %in% which(b_layers %in% lay_pre$idx:lay_ii$idx), ncol=s),
      1, all
    )
    # Do not count faces that are all made of pre-layer vertices, or
    #   all post-layer vertices.
    faces_btwn <- faces_btwn & (!(lay_pre$faces_complete)) & (!(lay_ii$faces_complete))
    faces2$rm[faces_btwn] <- TRUE

    # --------------------------------------------------------------------------
    # Make new faces. ----------------------------------------------------------
    # Strategy:
    #   * Start with an arbitrary vertex from layer A.
    #   * Make a face between that vertex, the closest (first) in layer B, and
    #     the second in layer A next in radial order).
    #   * Then, make a face between the second in layer A, the first in layer B,
    #   * and the second in layer B. These two faces are a square with a dividing
    #   * diagonal.
    #   * Continue until looped around.
    #   * Add triangles evenly distributed to account for different number
    #   * of vertices.
    # --------------------------------------------------------------------------

    # Multiply adjacency matrices between all in-between layers to get
    #   pseudo-adjacency matrix between vertices in the post- and pre- layers.
    adj <- b_adjies[lay_pre$idx:lay_ii$idx]
    adj <- Reduce("%*%", adj)
    # ??? Not used yet.

    V_max <- max(lay_ii$V2, lay_pre$V2)
    v_ii <- lay_ii$verts[order(lay_ii$rad_order)][seq(1, lay_ii$V1, lay_k[ii+1])]
    v_ii <- v_ii[as.integer(ceiling(seq(1, V_max) * (lay_ii$V2 / V_max)))]
    v_pre <- lay_pre$verts[order(lay_pre$rad_order)][seq(1, lay_pre$V1, lay_k[ii])]
    v_pre <- v_pre[as.integer(ceiling(seq(1, V_max) * (lay_pre$V2 / V_max)))]

    one_back <- function(x){c(x[2:length(x)], x[1])}
    cost_a <- mean(apply(vertices[v_ii,] - vertices[one_back(v_pre),], 1, norm))
    cost_b <- mean(apply(vertices[one_back(v_ii),] - vertices[v_pre,], 1, norm))
    # Note: will remove degenerate faces after loop (v_ii/v_pre have repeats)
    if (cost_a < cost_b) {
      faces2$new <- rbind(faces2$new, cbind(v_ii, one_back(v_ii), v_pre))
      faces2$new <- rbind(faces2$new, cbind(one_back(v_ii), v_pre, one_back(v_pre)))
    } else {
      faces2$new <- rbind(faces2$new, cbind(v_pre, one_back(v_pre), v_ii))
      faces2$new <- rbind(faces2$new, cbind(one_back(v_pre), v_ii, one_back(v_ii)))
    }
  }

  # ----------------------------------------------------------------------------
  # Arrange the results. -------------------------------------------------------
  # ----------------------------------------------------------------------------

  # Remove degenerate faces.
  faces2$degen <- (faces2$new[,1] == faces2$new[,2])
  faces2$degen <- faces2$degen | (faces2$new[,1] == faces2$new[,3])
  faces2$degen <- faces2$degen | (faces2$new[,2] == faces2$new[,3])
  faces2$new <- faces2$new[!faces2$degen,]

  # Construct and return new mesh.
  #faces2$new <- rbind(faces[!faces2$rm,], faces2$new)
  verts_remaining <- unique(as.vector(faces2$new))
  verts2 <- rep(NA, V_)
  verts2[verts_remaining] <- 1:length(verts_remaining)
  faces2$new <- matrix(verts2[as.vector(faces2$new)], ncol=s)
  list(vertices=vertices[verts_remaining,], faces=faces2$new)
}


#' Apply Mask to Vertices and Faces
#'
#' Apply a binary mask to a set of vertices and faces.  Vertices not in the mask are removed,
#' and faces (triangles) with any vertices outside of the mask are removed.  Finally,
#' vertex numbering in the masked faces matrix is corrected.
#'
#' @inheritParams vertices_Param
#' @inheritParams faces_Param
#' @inheritParams mask_Param_vertices
#'
#' @return List containing masked vertices and faces matrices
#' 
#' @export
mask_vertices_faces <- function(vertices, faces, mask){

  # Number of vertices
  V <- nrow(vertices)

  # Check index of faces
  if(min(faces) == 0){
    faces <- faces + 1
  }

  mask <- as.numeric(mask)
  if(length(mask) != V | !is.vector(mask)){
    stop("Mask should be a vector of length V")
  }
  # Check only 0s and 1s
  values <- sort(unique(mask))
  if(! (min(values %in% 0:1)) ) stop("Mask should be composed of only 0s and 1s")

  inmask <- which(mask==1)

  # Apply mask to vertices
  vertices_new <- vertices[inmask,]

  ### Apply mask to faces (triangles)

  # Identify triangles where any vertex is outside of the mask
  faces <- faces[(faces[,1] %in% inmask) & (faces[,2] %in% inmask) & (faces[,3] %in% inmask),]

  # Re-number faces
  faces_new <- faces*0
  for(ii in 1:nrow(faces)){
    faces_new[ii,1] <- which(inmask == faces[ii,1])
    faces_new[ii,2] <- which(inmask == faces[ii,2])
    faces_new[ii,3] <- which(inmask == faces[ii,3])
  }

  # Return updated vertices and faces
  result <- list(vertices=vertices_new, faces=faces_new)
  return(result)
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
		v_sd <- sqrt(rowVars(BOLD, na.rm=TRUE))
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

#' Find nonzero elements in matrix
#'
#' Find nonzero element in a matrix using 2-means clustering
#'
#' @importFrom stats kmeans
#'
#' @param beta_est A vector or matrix of values from which values close to zero
#'  should be assigned a value of zero.
#'
#' @return A vector or matrix of the same dimension as beta_est in which values
#'  close to zero are assigned the value of zero. The closeness of a value to
#'  zero is found by performing two-means clustering on the absolute values of
#'  beta_est, and ...
#'
#' @export
#'
find_nonzero <- function(beta_est) {
  vector_beta <- c(beta_est)
  if(any(is.na(vector_beta))) vector_beta <- vector_beta[!is.na(vector_beta)]
  km_beta <- kmeans(abs(vector_beta),2)
  which_nonzero <- which.max(km_beta$centers[,1])
  keep_nonzero <- as.numeric(km_beta$cluster == which_nonzero)
  out <- beta_est
  out[!is.na(out)] <- out[!is.na(out)] * keep_nonzero
  return(out)
}

#' Boundary Mask
#'
#' Identify the vertices within `boundary_width` edges of the input mask. The
#'  faces must be triangular.
#'
#' @param faces a V x 3 matrix of integers. Each row defines a face by the index
#'  of three vertices.
#' @inheritParams mask_Param_vertices
#' @param boundary_width a positive integer. Vertices no more than this number
#'  of edges from any vertex in the input mask will be placed in the boundary mask.
#'
#' @return The boundary mask, a length-V logical vector. TRUE indicates vertices
#'  within the boundary mask.
#'
#' @keywords internal
boundary_mask <- function(faces, mask, boundary_width){
  s <- ncol(faces)
  v <- max(faces)
  # For quads, boundary_mask() would count opposite vertices on a face as
  #   adjacent--that's probably not desired.
  stopifnot(s == 3)

  stopifnot(boundary_width > 0)

  boundary_mask <- rep(FALSE, v)
  # Begin with the input mask.
  verts_adj_previous <- which(mask)
  for (ii in seq(1, boundary_width)) {
    # Identify vertices not in the mask, but adjacent to it.
    # Adjacency is defined by sharing a face.
    faces_nmask <- rowSums(matrix(faces %in% verts_adj_previous, ncol=s))
    faces_adj <- faces_nmask > 0 & faces_nmask < s
    verts_adj <- unique(as.vector(faces[faces_adj,]))
    verts_adj <- verts_adj[!(verts_adj %in% verts_adj_previous)]
    # Add those vertices to the boundary mask, and use them as the mask in
    #   the next iteration.
    boundary_mask[verts_adj] <- TRUE
    verts_adj_previous <- verts_adj
  }

  boundary_mask
}

#' Match user inputs to expected values
#'
#' Match each user input to an expected/allowed value. Raise a warning if either
#'  several user inputs match the same expected value, or at least one could not
#'  be matched to any expected value. \code{ciftiTools} uses this function to
#'  match keyword arguments for a function call. Another use is to match
#'  brainstructure labels ("left", "right", or "subcortical").
#'
#' @param user Character vector of user input. These will be matched to
#'  \code{expected} using \code{match.arg()}.
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
    }
  )

  invisible(NULL)
}


#' Create a mask based on vertices that are invalid
#'
#' @param data A list of sessions, where each session is a list with elements
#'  BOLD, design and nuisance. See \code{?create.session} and \code{?is.session}
#'  for more details.
#' @param meanTol,varTol Tolerance for mean and variance of each data location. Locations which
#'  do not meet these thresholds are masked out of the analysis. Defaults: \code{1e-6}.
#' @param verbose Print messages counting how many locations are removed?
#'
#' @importFrom matrixStats colVars
#' @return A logical vector indicating valid vertices
#'
#' @export
make_mask <- function(data, meanTol=1e-6, varTol=1e-6, verbose=TRUE){

  # For each BOLD data matrix,
  mask_na <- mask_mean <- mask_var <- rep(TRUE, ncol(data[[1]]$BOLD))
  for (ss in seq(length(data))) {
    dss <- data[[ss]]$BOLD
    # Mark columns with any NA or NaN values for removal.
    dss_na <- is.na(dss) | is.nan(dss)
    mask_na[apply(dss_na, 2, any)] <- FALSE
    # Calculate means and variances of columns, except those with any NA or NaN.
    # Mark columns with mean/var falling under the thresholds for removal.
    mask_mean[mask_na][colMeans(dss[,mask_na,drop=FALSE]) < meanTol] <- FALSE
    mask_var[mask_na][matrixStats::colVars(dss[,mask_na,drop=FALSE]) < varTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- if (any(!mask_na)) { "additional locations" } else { "locations" }
    warn_part2 <- if (length(data) > 1) { " in at least one scan.\n" } else { ".\n" }
    if (any(!mask_na)) {
      cat("\t", sum(!mask_na), paste0("locations removed due to NA or NaN values", warn_part2))
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat("\t", sum(!mask_mean2), warn_part1, paste0("removed due to low mean", warn_part2))
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat("\t", sum(!mask_var2), warn_part1, paste0("removed due to low variance", warn_part2))
    }
  }

  # Return composite mask.
  mask_na & mask_mean & mask_var
}

#' Transform vector data to an image
#'
#' This fills in parts of a template with values from \code{vec_data}.
#'
#' @param vec_data A V by p matrix, where V is the number of voxels within a
#'   mask and p is the number of vectors to transform into matrix images
#' @param template_image A binary matrix in which V entries are 1 and the rest
#'   of the entries are zero
#'
#' @return A list of masked values from \code{vec_data}
#' 
#' @export
vec2image <- function(vec_data, template_image) {
  each_col <- sapply(split(vec_data, col(vec_data)), function(vd) {
    out <- template_image
    out[out == 1] <- vd
    out[out == 0] <- NA
    return(out)
  }, simplify = F)
  return(each_col)
}
