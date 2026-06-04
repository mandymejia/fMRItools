#' Find best matches for each network based on Dice overlap
#'
#' For each network (map or parcel) in A, it finds the network in B with which
#'  it has the highest Dice coefficient. (Some networks in B will have zero
#'  matches, and some will have more than one.)
#'
#' If the resolutions of A and B differ, it will try to resample B to match
#'  the resolution of A.
#'
#' @param A The networks to find matches for. Can be a CIFTI file path, a
#'  \code{"xifti"} object, a numeric matrix (locations by networks, and can be
#'  continuous or binary), or a locations-length vector giving the network index
#'  each location belongs to. 
#' @param B The networks to match to. Can be a CIFTI file path, a
#'  \code{"xifti"} object, a numeric matrix (locations by networks, and can be
#'  continuous or binary), a locations-length vector giving the network index
#'  each location belongs to, or the name of a parcellation included in 
#'  \code{ciftiTools}. Default: \code{"Yeo_7"}. 
#' @param B_zero_drop If B is a parcellation (i.e. a vector of network indices),
#'  and the lowest parcel index is 0, drop that level? Sometimes the zero-valued
#'  parcel is used for the medial wall or other invalid locations. Default: 
#'  \code{TRUE}. If \code{FALSE}, treat it as a candidate for matching just like
#'  all the other parcels.
#' @param brainstructures The brainstructures to use, if \code{A} and/or
#'  \code{B} are CIFTI data. Default: \code{"existing"} to use the
#'  brainstructures present in \code{A} (or in \code{B} if \code{A} is a matrix).
#' 
#' @return A numeric vector giving the indices of \code{B} networks that are the best
#'  matches for each network in \code{A}. \code{order(match_nets(A, B))} would give an
#'  indexing of \code{A} that rearranges its networks according to the order of
#'  networks in \code{B}.
#'
#' @export
#'
match_nets <- function(A, B="Yeo_7", B_zero_drop=TRUE, brainstructures="existing") {

  # If either A or B are not logical/numeric vectors/matrices
  if ((is.list(A) || is.character(A)) || (is.list(B) || is.character(B))) {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it", call. = FALSE)
    }

    # Check brainstructures match, if `A` and `B` are CIFTI.
    bs_all <- c("cortex_left", "cortex_right", "subcort")
    if (is_1(A, "character")) {
      A <- ciftiTools::read_cifti(A, brainstructures=brainstructures)
    }
    if (ciftiTools::is.xifti(A, message=FALSE)) {
      bs_here <- bs_all[vapply(bs_all, function(q){!is.null(A$data[[q]])}, FALSE)]
    } else {
      bs_here <- "existing" # determine based on B
    }

    # Load/read B
    if (is_1(B, "character")) {
      if (B %in% c("Schaefer_100", "Schaefer_400", "Schaefer_1000", "Yeo_7", "Yeo_17")) {
        B <- ciftiTools::load_parc(B)
      } else {
        B <- ciftiTools::read_cifti(B, brainstructures = bs_here)
      }
    }

    # Try to handle mismatch of number of rows. 
    # Several reasons: brainstructures, medial wall, resolution.
    if (nrow(A) != nrow(B)) {
      # Set medial wall values to zero. Handles discrepancies in medial wall shape/size.
      # [NOTE] Assumes that if `A` is a matrix, it includes the medial wall.
      if (ciftiTools::is.xifti(A, message=FALSE)) { A <- ciftiTools::move_from_mwall(A, 0) }
      if (ciftiTools::is.xifti(B, message=FALSE)) { B <- ciftiTools::move_from_mwall(B, 0) }

      # `bs_here`: what's in A, or B, if A is not a xifti.
      # `brainstructures`: what's requested, confirmed in A if xifti, and wanted in B.
      if (ciftiTools::is.xifti(B, message=FALSE)) {
        if (bs_here == "existing") {
          bs_here <- bs_all[vapply(bs_all, function(q){!is.null(B$data[[q]])}, FALSE)]
        }
        if (is_1(brainstructures, "character") && brainstructures=="existing") {
          message("Assuming `A` has these brainstructures: ", paste(bs_here, collapse=", "), ".")
          brainstructures <- bs_here
        }
        bs_missing <- setdiff(brainstructures, bs_here)
        if (length(bs_missing) > 0) {
          stop("These brainstructures are missing from `B`: ",
               paste(bs_missing, collapse=", "), "."
          )
        }

        # Remove bs in B not in A/brainstructures. If any bs in A/brainstructures not in B, raise error.
        bs_B <- bs_all[vapply(bs_all, function(q){!is.null(B$data[[q]])}, FALSE)]
        if (!(length(brainstructures) == length(bs_B)) || !all(order(bs_B) == order(brainstructures))) {
          B <- ciftiTools::remove_xifti(B, remove=setdiff(bs_B, brainstructures))
          bs_B <- bs_all[vapply(bs_all, function(q){!is.null(B$data[[q]])}, FALSE)]
          bs_missing <- setdiff(brainstructures, bs_B)
          if (length(bs_missing) > 0) {
            stop("These brainstructures in `A` are missing from `B`: ",
                 paste(bs_missing, collapse=", "), "."
            )
          }
        }

        # Resample if necessary and possible.
        if (nrow(A) != nrow(B) && ("cortex_left" %in% brainstructures || "cortex_right" %in% brainstructures)) {
          message("`A` and `B` have different numbers of locations (rows). Will try to resample the cortex.")
          A_cortex_res <- max(nrow(A$data$cortex_left), nrow(A$data$cortex_right), na.rm=TRUE)
          B_cortex_res <- max(nrow(B$data$cortex_left), nrow(B$data$cortex_right), na.rm=TRUE)
          message("`A` total number of rows: ", nrow(A))
          message("`B` total number of rows: ", nrow(B))
          message("`A` cortex vertices: ", A_cortex_res)
          message("`B` cortex vertices: ", B_cortex_res)
          B <- ciftiTools::resample_xifti(B, resamp_res=A_cortex_res)
          message("new `B` number of rows: ", nrow(B))
        }
      }
    }
  }

  A <- as.matrix(A)
  B <- as.matrix(B)

  # Handle `NA` values in matrix case.
  if (ncol(A) > 1 && any(is.na(A))) {
    na_sums <- rowSums(is.na(A))
    if (!all(na_sums %in% c(0, ncol(A)))) {
      stop("Not allowed: locations `NA` for only some columns of `A`.")
    }
    A[na_sums == ncol(A),] <- 0
  }
  if (ncol(B) > 1 && any(is.na(B))) {
    na_sums <- rowSums(is.na(B))
    if (!all(na_sums %in% c(0, ncol(B)))) {
      stop("Not allowed: locations `NA` for only some columns of `B`.")
    }
    B[na_sums == ncol(B),] <- 0
  }

  # One-hot encode, and handle `NA` values, in vector case.
  ## A.
  if (ncol(A) == 1 && all(c(A) == round(c(A)), na.rm=TRUE)) {
    first_idx <- min(A, na.rm=TRUE)
    first_drop <- FALSE

    na_mask <- is.na(A)
    if (any(na_mask)) {
      first_drop <- TRUE
      A[na_mask] <- first_idx - 1 # put `NA` in the drop level.
    }

    minA <- min(A)
    if (minA < 1) {
      # Warn users if A indices will shift, but don't discuss the temp drop level for NAs.
      if (!(minA==0 && first_drop==TRUE)) {
        message("Shifting indices of A by adding ", -minA+1-first_drop, " to each so the first index is 1.")
      }
      A <- as.numeric(A) -minA+1
      idx2 <- sort(unique(A))
      if (first_drop) { idx2 <- idx2[-1] }
      if (!(minA==0 && first_drop==TRUE)) {
        if (any(diff(idx2) > 1)) {
          message("The new indices are: ", paste(idx2, collapse=", "))
        }
      }
    }
    A <- diag(max(A))[A, ] > 0 # the one-hot matrix. yes, this works!
    if (first_drop) {
      A <- A[,-1]
    }
  }

  ## B.
  ## Different than A in that any index shifting is handled silently.
  temp_shift <- 0
  if (ncol(B) == 1 && all(c(B) == round(c(B)), na.rm=TRUE)) {
    first_idx <- min(B, na.rm=TRUE)
    first_drop <- FALSE

    if (first_idx == 0 && B_zero_drop) {
      first_drop <- TRUE
      first_idx <- 1
    }

    na_mask <- is.na(B)
    if (any(na_mask)) {
      first_drop <- TRUE
      B[na_mask] <- first_idx - 1 # put `NA` in the drop level.
    }

    minB <- min(B)
    if (minB < 1) {
      temp_shift <- -minB+1 # quietly shift, but will return indices w/ original values
      B <- as.numeric(B) + temp_shift
    }
    B <- diag(max(B))[B, ] > 0 # the one-hot matrix. yes, this works!
    if (first_drop) {
      B <- B[,-1]
      temp_shift <- temp_shift - 1
    }
  }

  if (nrow(A) != nrow(B)) { stop("`A` has ", nrow(A), " rows but `B` has ", nrow(B), " rows.") }
  dice_mat <- dice_coef(A, B)
  best_match <- apply(dice_mat, 1, which.max)
  best_match - temp_shift
}