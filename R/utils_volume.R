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
