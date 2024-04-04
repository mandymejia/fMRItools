#' Un-vectorize correlation matrix
#'
#' Converted a vectorized lower triangular correlation matrix back to its full
#'  matrix form.
#'
#' @param x_diag The vectorized lower triangular correlation matrix.
#' @param diag_val The value to put on the diagonals of the correlation matrix.
#'  Default: \code{NA}.
#' @param names (Optional) row/column names.
#' @param newOrder (Optional) new index order.
#' @param lowerOnly Set the upper triangle to \code{NA}? Default: \code{FALSE}.
#'
#' @return A square correlation matrix
#' @keywords internal
cor_mat <- function(x_diag, diag_val=NA, names=NULL, newOrder=NULL, lowerOnly=FALSE) {
  d <- 1/2 + sqrt(1/4 + 2*length(x_diag))
  mat <- diag(d) * diag_val
  mat[upper.tri(mat)] <- x_diag
  mat <- t(mat)
  mat[upper.tri(mat)] <- x_diag
  if (!is.null(names)) {
    stopifnot(length(names)==d); rownames(mat) <- colnames(mat) <- names
  }
  if (!is.null(newOrder)) {
    stopifnot(all(sort(newOrder) == seq(d))); mat <- mat[newOrder,rev(newOrder)]
  }
  if (lowerOnly) { mat[seq(d), rev(seq(d))][lower.tri(mat)] <- NA }
  mat
}

#' Color palette
#'
#' Color palettes for fMRI data analysis tasks
#'
#' @param pal \code{"div_Beach"} (default; blue to white to red), 
#'  \code{"seq_Beach"} (white to red), or 
#'  \code{"seq_Beach2"} (blue to white).
#' @return A data.frame with two columns: \code{"col"} has the hex code of color,
#' and \code{"val"} has the placement of colors between zero and one.
#' @export
color_palette <- function(pal="div_Beach") {
  switch(pal,
    div_Beach = ciftiTools::expand_color_pal(data.frame(
      color = c(
        "#1a188f", "#5e5eff", "#78bbff", "#9bf5ff",
        "#e1f7e9", "#fbfff9", "#f5f5da",
        "#fafa89", "#ffa500", "#ff2424", "#680000"
      ),
      value = c(0, .1, .2, .32, .42, .5, .58, .68, .8, .9, 1)
    ), COLOR_RES=400)$color,
    seq_Beach = ciftiTools::expand_color_pal(data.frame(
      color = c(
        "#fbfff9", "#f5f5da",
        "#fafa89", "#ffa500", "#ff2424", "#680000"
      ),
      value = c(.5, .58, .68, .8, .9, 1) * 2 - 1
    ), COLOR_RES=400)$color,
    seq_Beach2 = ciftiTools::expand_color_pal(data.frame(
      color = rev(c(
        "#1a188f", "#5e5eff", "#78bbff", "#9bf5ff",
        "#e1f7e9", "#fbfff9"
      )),
      value = rev(c(0, .1, .2, .32, .42, .5) * 2)
    ), COLOR_RES=400)$color
  )
}

#' image.scale
#'
#' image.scale. Source: r-bloggers.com/2013/12/new-version-of-image-scale-function/
#'
#' @param z,zlim,col,breaks,axis.pos,add.axis,... The arguments.
#' @return Plots the image scale.
#' @keywords internal
image.scale <- function(z, zlim, col = color_palette("div_Beach"),
  breaks, axis.pos=1, add.axis=TRUE, ...){

  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed. Please install it", call. = FALSE)
  }

  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }

  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }

  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }

  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }

  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}

  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)

  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      graphics::polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }

    if(axis.pos %in% c(2,4)){
      graphics::polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }

 graphics::box()

 if(add.axis) {graphics::axis(axis.pos)}
}

#' Plot FC
#'
#' Plot a functional connectivity matrix.
#'
#' @param FC The functional connectivity matrix, a square numeric matrix with
#'  values between -1 and 1.
#' @param zlim The minimum and maximum range of the color scale. Default:
#'  \code{c(-1, 1)}. If in descending order, the color scale will be reversed.
#' @param diag_val Set to \code{NA} for white, \code{1}, or \code{NULL} 
#'  (default) to not modify the diagonal values in \code{FC}.
#' @param title (Optional) Plot title.
#' @param cols Character vector of colors for the color scale. Default:
#'  \code{color_palette("div_Beach")}.
#' @param cleg_ticks_by Spacing between ticks on the color legend. Default:
#'  \code{0.5}.
#' @param cleg_digits How many decimal digits for the color legend. Default:
#'  \code{1}.
#' @param labels A character vector of length \code{length(lines)+1}, giving
#'  row/column labels for the submatrices divided by \code{lines}. If 
#'  \code{NULL} (default), do not add labels.
#' @param lines Add lines to divide the FC matrix into submatrices? Default:
#'  delineate each individual row and column. 
#' @param lines_col,lines_lwd Color and line width of the \code{lines}. Default:
#'  black lines of width \code{1}.
#' @param cex Text size. Default: \code{0.8}.
#'
#' @export
plot_FC <- function(
  FC, zlim=c(-1,1),
  diag_val=NULL,
  title="FC matrix",
  cols=color_palette("div_Beach"),
  cleg_ticks_by=0.5, cleg_digits=1,
  labels = NULL,
  lines = seq(nN),
  lines_col = 'black',
  lines_lwd = 1,
  cex = 0.8
  ){

  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package \"grDevices\" needed. Please install it", call. = FALSE)
  }
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed. Please install it", call. = FALSE)
  }

  old_par <- graphics::par(no.readonly = TRUE)

  # Prep FC matrix -----
  nN <- ncol(FC)
  if (!is.null(diag_val)) { diag(FC) <- diag_val }
  # Truncate values to within zlim.
  FC[] <- pmax(FC[], zlim[1])
  FC[] <- pmin(FC[], zlim[2])

  # Prep colors -----
  # Reverse color scale if zlim is in descending order.
  if (zlim[2] < zlim[1]) {
    cols <- rev(cols)
    zlim <- rev(zlim)
  }
  # Make the color scale higher-resolution, if necessary.
  color_res <- 11 #401
  if (length(cols) < color_res) {
    cols <- grDevices::colorRampPalette(cols, space="Lab")(color_res)
  }
  cleg_ticks <- format(round(
    seq(zlim[1], zlim[2], cleg_ticks_by), cleg_digits
  ), scientific=FALSE, digits = cleg_digits, nsmall = cleg_digits)

  # Plot -----
  graphics::layout(matrix(c(1,2,0,3), nrow=2, ncol=2), widths=c(5,1.2), heights=c(1.2,5))

  ### Title -----
  graphics::par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  graphics::text(x = 0.5, y = 0.3, title, cex = 2, col = "black")

  ### FC plot -----
  graphics::par(mar=c(1,2,0,1))
  graphics::image(
    seq(nN), seq(nN), t(FC[rev(seq(nN)),]), col=cols,
    zlim=zlim,
    xaxt="n", yaxt="n", ylab="", xlab=""
  )
  ##### Lines -----
  graphics::abline(h=nN-lines+0.5, v=lines+0.5, col = lines_col, lwd=lines_lwd)
  ##### Labels -----
  if(!is.null(labels)){
    nK <- length(lines) + 1
    if(length(labels) == nK){
      for (k in seq(nK)) {
        k_start <- if(k == 1) { 0 } else { lines[k-1] }
        k_end <- if(k == nK) { nN } else { lines[k] }
        k_at <- k_start + (k_end-k_start)/2 + 0.5
        graphics::mtext(labels[k], side = 3, at = k_at, line=0, font=2, cex=cex)
        graphics::mtext(labels[k], side = 2, at = nN - k_at + 1, line=0, font=2, cex=cex)
      }
    } else {
      warning("Ignoring `labels`. `lines` divides the FC matrix into ", nK, " rows/cols, which should match the number of labels.")
    }
  }

  ### Color scale -----
  graphics::par(mar=c(1, 0.7, 0, 4))
  image.scale(
    FC, col=cols,
    zlim=zlim,
    axis.pos=4, add.axis=FALSE
  )
  graphics::axis(4, at=cleg_ticks, las=2, labels=cleg_ticks)
  graphics::abline(h=cleg_ticks)

  graphics::par(old_par)
}