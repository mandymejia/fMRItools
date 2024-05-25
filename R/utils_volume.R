#' Pad 3D Array
#'
#' Pad a 3D array by a certain amount in each direction, along each dimension.
#'  This operation is like the opposite of cropping.
#'
#' @param x A 3D array, e.g.
#'  \code{unvec_vol(xifti$data$subcort, xifti$meta$subcort$mask)}.
#' @param padding A \eqn{3 \times 2} matrix indicating the number of
#'  slices to add at the beginning (first column) and end (second column) of
#'  each of dimension, e.g. \code{xifti$meta$subcort$mask_padding}.
#' @param fill Value to pad with. Default: \code{NA}.
#'
#' @return The padded array
#'
#' @export
#'
#' @examples
#' x <- array(seq(24), dim=c(2,3,4))
#' y <- pad_vol(x, array(1, dim=c(3,2)), 0)
#' stopifnot(all(dim(y) == dim(x)+2))
#' stopifnot(sum(y) == sum(x))
#' z <- crop_vol(y)$data
#' stopifnot(identical(dim(x), dim(z)))
#' stopifnot(max(abs(z - x))==0)
pad_vol <- function(x, padding, fill=NA){
  # Argument checks.
  stopifnot(length(dim(x)) == 3)
  stopifnot(length(dim(padding)) == 2)
  stopifnot(all(dim(padding) == c(3,2)))
  class(padding) <- "numeric"
  stopifnot(all(vapply(c(padding), is_integer, nneg=TRUE, TRUE)))
  stopifnot(length(fill) == 1)

  # Get dimensions of new array.
  new_dim <- dim(x) + padding[,1] + padding[,2]

  # Get new array.
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

#' Sum of each voxel's neighbors
#'
#' For each voxel in a 3D logical or numeric array, sum the values of the six
#'  neighboring voxels.
#'
#' Diagonal voxels are not considered adjacent, i.e. the voxel at (0,0,0) is not
#'  adjacent to the voxels at (1,1,0) or (1,1,1), although it is adjacent to
#'  (1,0,0).
#'
#' @param arr The 3D array.
#' @param pad In order to compute the sum, the array is temporarily padded along
#'  each edge with the value of \code{pad}. \code{0} (default) will mean that
#'  edge voxels reflect the sum of 3-5 neighbors whereas non-edge voxels reflect
#'  the sum of 6 neighbors. An alternative is to use a value of \code{NA} so
#'  that edge voxels are \code{NA}-valued because they did not have a complete
#'  set of six neighbors. Perhaps another option is to use \code{mean(arr)}.
#'
#' @return An array with the same dimensions as \code{arr}. Each voxel value
#'  will be the sum across the immediate neighbors. If \code{arr} was a logical
#'  array, this value will be between 0 and 6.
#'
#' @export
sum_neighbors_vol <- function(arr, pad=0){
  # Argument checks.
  stopifnot(length(dim(arr)) == 3)
  class(arr) <- "numeric"
  stopifnot(length(pad) == 1)

  # Pad array by one voxel at the start and end of each dim.
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

#' Erode 3D mask
#'
#' Erode a volumetric mask by a certain number of voxel layers. For each layer,
#'  any in-mask voxel adjacent to at least one out-of-mask voxel is removed
#'  from the mask.
#'
#' Diagonal voxels are not considered adjacent, i.e. the voxel at (0,0,0) is not
#'  adjacent to the voxels at (1,1,0) or (1,1,1), although it is adjacent to
#'  (1,0,0).
#'
#' @param vol The 3D array to erode. The mask to erode is defined by all values
#'  not in \code{out_of_mask_val}.
#' @param n_erosion The number of layers to erode the mask by. Default:
#'  \code{1}.
#' @param out_of_mask_val A voxel is not included in the mask if and only if its
#'  value is in this vector. The first value of this vector will be used to
#'  replace eroded voxels. Default: \code{NA}. If \code{vol} is simply a logical
#'  array with \code{TRUE} values for in-mask voxels, use
#'  \code{out_of_mask_val=FALSE}.
#'
#' @return The eroded \code{vol}. It is the same as \code{vol}, but eroded
#'  voxels are replaced with \code{out_of_mask_val[1]}.
#'
#' @export
erode_mask_vol <- function(vol, n_erosion=1, out_of_mask_val=NA){
  # Argument checks.
  stopifnot(length(dim(vol)) == 3)
  stopifnot(is_integer(n_erosion, nneg=TRUE))
  if (n_erosion == 0) { return(vol) }

  # Erode.
  mask <- !array(vol %in% out_of_mask_val, dim=dim(vol))
  for (ii in seq(n_erosion)) {
    to_erode <- mask & (sum_neighbors_vol(mask, pad=1) < 6)
    mask[to_erode] <- FALSE
    vol[to_erode] <- out_of_mask_val[1]
  }
  vol
}

#' Dilate 3D mask
#'
#' Dilate a volumetric mask by a certain number of voxel layers. For each layer,
#'  any out-of-mask voxel adjacent to at least one in-mask voxel is added to the
#'  mask.
#'
#' Diagonal voxels are not considered adjacent, i.e. the voxel at (0,0,0) is not
#'  adjacent to the voxels at (1,1,0) or (1,1,1), although it is adjacent to
#'  (1,0,0).
#'
#' @param vol The 3D array to dilate. The mask to dilate is defined by all
#'  values not in \code{out_of_mask_val}.
#' @param n_dilate The number of layers to dilate the mask by. Default:
#'  \code{1}.
#' @param out_of_mask_val A voxel is not included in the mask if and only if its
#'  value is in this vector. Default: \code{NA}. If \code{vol} is simply a
#'  logical array with \code{TRUE} values for in-mask voxels, use
#'  \code{out_of_mask_val=FALSE}.
#' @param new_val Value for voxels newly added to the mask. Default: \code{1}.
#'  If \code{vol} is simply a logical array with \code{TRUE} values for
#'  in-mask voxels, use \code{new_val=1}.
#'
#' @return The dilated \code{vol}. It is the same as \code{vol}, but dilated
#'  voxels are replaced with \code{new_val}.
#'
#' @export
dilate_mask_vol <- function(vol, n_dilate=1, out_of_mask_val=NA, new_val=1){
  # Argument checks.
  stopifnot(length(dim(vol)) == 3)
  stopifnot(is_integer(n_dilate, nneg=TRUE))
  stopifnot(is_1(new_val, "numeric") | is_1(new_val, "logical"))
  if (n_dilate == 0) { return(vol) }

  # Dilate.
  mask <- !array(vol %in% out_of_mask_val, dim=dim(vol))
  for (ii in seq(n_dilate)) {
    to_dilate <- (!mask) & (sum_neighbors_vol(mask, pad=0) > 0)
    mask[to_dilate] <- TRUE
    vol[to_dilate] <- new_val
  }
  vol
}

#' Convert coordinate list to 3D array
#'
#' Converts a sparse coordinate list to its non-sparse volumetric representation.
#'
#' @param coords The sparse coordinate list. Should be a \code{"data.frame"} or
#'  matrix with voxels along the rows and three or four columns. The first three
#'  columns should be integers indicating the spatial coordinates of the voxel.
#'  If the fourth column is present, it will be the value used for that voxel.
#'  If it is absent, the value will be \code{TRUE} or \code{1} if \code{fill}
#'  is not one of those values, and \code{FALSE} or \code{0} if \code{fill} is.
#'  The data type will be the same as that of \code{fill}. The fourth column
#'  must be logical or numeric.
#' @param fill Logical or numeric fill value for the volume. Default:
#'  \code{FALSE}.
#'
#' @return The volumetric data
#'
#' @export
#'
coordlist_to_vol <- function(coords, fill=FALSE){
  # Check arguments.
  stopifnot(length(fill) == 1)
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
#' @param x The numeric 3D array to crop.
#' @return A list of length two: \code{"data"}, the cropped array, and
#'  \code{"padding"}, the number of slices removed from each edge of each
#'  dimension.
#'
#' @export
#'
crop_vol <- function(x) {
  # Argument checks.
  d <- length(dim(x))
  stopifnot(d == 3)
  if (all(unique(as.vector(x)) %in% c(NA, 0))) { stop("The array is empty.") }

  # For each dimension, get the first and last non-empty slice.
  # This chunk of code is generalizable to arrays with up to 15 dimensions.
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

  # Crop.
  # This can easily be generalized to support arrays with number of dimensions
  #   other than 3. But for now, it only works for 3D arrays.
  x <- x[!empty_slice[[1]], !empty_slice[[2]], !empty_slice[[3]]]

  list(data=x, padding=padding)
}

#' Get coordinates of each voxel in a mask
#'
#' Made for obtaining voxel locations in 3D space from the subcortical metadata
#'  of CIFTI data: the volumetric mask, the transformation matrix and the
#'  spatial units.
#'
#' @param mask 3D logical mask
#' @param trans_mat Transformation matrix from array indices to spatial
#'  coordinates.
#' @param trans_units Units for the spatial coordinates (optional).
#' @return A list: \code{coords} and \code{trans_units}.
#'
#' @export
#'
vox_locations <- function(mask, trans_mat, trans_units=NULL){
  # Argument checks.
  stopifnot(dim(mask) == 3)
  stopifnot(class(mask) == "logical")
  stopifnot(dim(class(trans_mat)) == 2)
  stopifnot(class(trans_mat) == "numeric")
  stopifnot(dim(trans_mat)[1] == 4)

  coords <- (cbind(which(mask, arr.ind=TRUE), 1) %*% trans_mat)[,seq(3)]
  list(coords = coords, trans_units = trans_units)
}

#' Convert vectorized data back to volume
#'
#' Un-applies a mask to vectorized data to yield its volumetric representation.
#'  The mask and data should have compatible dimensions: the number of rows in
#'  \code{dat} should equal the number of locations within the \code{mask}.
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
unvec_vol <- function(dat, mask, fill=NA) {

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
