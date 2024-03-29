% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_volume.R
\name{unvec_vol}
\alias{unvec_vol}
\title{Convert vectorized data back to volume}
\usage{
unvec_vol(dat, mask, fill = NA)
}
\arguments{
\item{dat}{Data matrix with locations along the rows and measurements along
the columns. If only one set of measurements were made, this may be a
vector.}

\item{mask}{Volumetric binary mask. \code{TRUE} indicates voxels inside the
mask.}

\item{fill}{The value for locations outside the mask. Default: \code{NA}.}
}
\value{
The 3D or 4D unflattened volume array
}
\description{
Un-applies a mask to vectorized data to yield its volumetric representation.
The mask and data should have compatible dimensions: the number of rows in
\code{dat} should equal the number of locations within the \code{mask}.
}
