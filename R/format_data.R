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
          ROI_noise[[ii]] <- erode_mask_vol(ROI_noise[[ii]], noise_erosion[[ii]], c(-1, 0, NA, NaN))
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
