#' Plot FC with ggplot2
#'
#' Plot a FC matrix.
#'
#' @param mat The FC matrix
#' @param colFUN A \code{scale_fill} function. If \code{NULL}, use a red-blue
#'  diverging palette if there are negative values in \code{mat}, and a
#'  white-blue sequential palette otherwise.
#' @param title The plot title. Default: \code{"FC Matrix"}.
#' @param legTitle The legend title. Default: \code{"FC"}.
#' @param group_divs,group_cols Split the FC matrix into groups of contiguous
#'  rows/columns? Use \code{group_divs} to indicate the index of the first
#'  element in each group. For example, if the groups are 1-8, 9-15, and 16 to
#'  25, set \code{group_divs=c(1, 9, 16)}. Groups will be indicated by a color
#'  bar along the left and top sides of the plot. Use \code{group_cols} to
#'  change the colors. By default, \code{group_divs} is \code{NULL} (no groups).
#' @param labs Labels to place along the y-axis? If \code{is.null(group_divs)},
#'  this should be a character vector labeling each row of \code{mat}, or a
#'  dot dot dot. If \code{!is.null(group_divs)}, this may alternatively be a
#'  character vector labeling each group.
#' @param uppertri_means Use the upper triangle of the FC matrix to show the
#'  mean of each group indicated by \code{group_divs}? Default: \code{TRUE}.
#'  Has no effect if \code{group_divs} is \code{NULL}.
#' @param divColor Color of group dividers and outline of the FC matrix.
#'  Default: \code{"black"}.
#' @param lim Length-two numeric vector indicating the limits for the color bar.
#'  Values beyond \code{lim} will be clipped. If \code{NULL} (default), the
#'  limits will be set to the largest magnitude value in \code{mat} (no
#'  clipping).
#' @param diagVal On-diagonal values will be set to this value.
#'  (\code{uppertri_means} are calculated before \code{diagVal} is used.)
#' @param labs_margin_y,labs_margin_x Margin value for labels. Default:
#'  \code{pmin(0, ncol(mat)/10-10)} for the y-axis and \code{0} for the x-axis.
#' @return The plot
#' @export
plot_FC_gg <- function(
  mat,
  colFUN=NULL,
  title="FC Matrix", legTitle="FC",
  group_divs=NULL, group_cols=RColorBrewer::brewer.pal(8, "Set2"),
  labs=NULL, uppertri_means=TRUE, divColor="black", lim=NULL, diagVal=1,
  labs_margin_y=pmin(0, ncol(mat)/10-10), labs_margin_x=0
  ){

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package \"RColorBrewer\" needed for FC plot. Please install it", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" needed for FC plot. Please install it", call. = FALSE)
  }
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
    stop("Package \"ggcorrplot\" needed for FC plot. Please install it", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for FC plot. Please install it", call. = FALSE)
  }

  nD <- nrow(mat)
  stopifnot(nD==ncol(mat))

  dw <- .004 * nD

  use_groups <- !is.null(group_divs)
  if (use_groups) {
    stopifnot(all(diff(group_divs) > 0))
    stopifnot(group_divs[1] > 0)
    if (group_divs[1] != 1) { group_divs <- c(1, group_divs) }
    stopifnot(group_divs[length(group_divs)] < nD+2)
    if (group_divs[length(group_divs)] != nD+1) { group_divs <- c(group_divs, nD+1) }
    group_divs2 <- group_divs[seq(2, length(group_divs)-1)]
    gg_pdv <- data.frame(xmin=group_divs2-.5-(dw/2), xmax=group_divs2-.5+(dw/2))
    gg_pdv$Var1 <- gg_pdv$Var2 <- 0
    gdiv_cols <- as.numeric(cut(seq(nD)+1, c(-Inf, group_divs2, Inf)))
  }

  # Take network-network means in upper tri
  if (uppertri_means && use_groups) {
    mat2 <- mat
    for (jj in seq(length(group_divs)-1)) {
      for (kk in seq(length(group_divs)-1)) {
        mat2[seq(group_divs[jj], group_divs[jj+1]-1), nD+1 - seq(group_divs[kk], group_divs[kk+1]-1)] <- mean(
          mat2[seq(group_divs[jj], group_divs[jj+1]-1), nD+1 - seq(group_divs[kk], group_divs[kk+1]-1)], na.rm=TRUE
        )
      }
    }
    mat[seq(nD), rev(seq(nD))][lower.tri(mat)] <- mat2[seq(nD), rev(seq(nD))][lower.tri(mat)]
  }

  if (!is.null(diagVal)) { diag(mat) <- diagVal }

  # Limits
  if (!is.null(lim)) {
    if (length(lim)==1) { lim <- c(if(min(mat[])>=0) {0} else {-lim}, lim) }
    stopifnot(length(lim)==2)
  } else {
    lim <- max(abs(c(mat[upper.tri(mat)], mat[lower.tri(mat)])))
    lim <- c(if(min(mat[])>=0) {0} else {-lim}, lim)
  }
  #mat[] <- pmax(lim[1], pmin(lim[2], mat)) # scales::squish instead

  if (is.null(colFUN)) {
    use_seq <- lim[1] >= 0
    gvals <- c(
      "#5d0928", "#892041", "#b63b59", "#d26767", "#eb8e74", "#f4b794", "#fcdcba",
      "#fffee6",
      "#c7e7d8", "#99cbcf", "#68aec6", "#478db7", "#226ca7", "#1b4984", "#132560"
    )
    if (use_seq) { gvals <- gvals[seq(8, 15)] }
    colFUN <- function(limits=c(if(use_seq) {0} else {-1},1), ...) {
      ggplot2::scale_fill_gradientn(colours=gvals, limits=limits, ...)
    }
  }

  lim_expand <- ifelse(use_groups, nD*.1, 0)

  # Handle labs
  if (is.null(labs)) {
    breaks <- NULL
  } else if (use_groups && length(labs)==length(group_divs)-1) {
    breaks <- group_divs[seq(length(group_divs)-1)]
  } else if (length(labs)==nD) {
    breaks <- seq(nD)
  } else {
    stop("`labs` length should match the number of rows or groups in the FC matrix.")
  }

  mat <- mat[,rev(seq(ncol(mat)))] # couldn't see how to rev y-axis w/ coord_equal also happening

  plt <- ggcorrplot::ggcorrplot(mat, outline.color = "#00000000", title=title, digits=12)
  plt$scales$scales <- list() # stops warning about replacing scales
  plt$coordinates$default <- TRUE # stops warning about replacing coordinates

  plt <- plt +
    ggplot2::scale_x_continuous(breaks=nD-breaks+1, labels=rev(labs)) +
    ggplot2::scale_y_continuous(breaks=rev(nD-breaks+1), labels=rev(labs)) +
    ggplot2::coord_equal(
      xlim=.5+c(-1-lim_expand-dw, nD+dw),
      ylim=.5+c(-dw, nD+1+lim_expand+dw),
      expand=FALSE) +
    colFUN(
      c(lim[1], lim[2]), guide=ggplot2::guide_colorbar(ticks.colour = divColor, ticks.linewidth = 1),
      labels = function(x){gsub("0.", ".", x, fixed=TRUE)}, na.value="red", oob = scales::squish
    ) +
    ggplot2::labs(fill=legTitle) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(margin=ggplot2::margin(r=labs_margin_y)),
      axis.text.x = ggplot2::element_text(angle=45, margin=ggplot2::margin(t=labs_margin_x)),
      #axis.ticks.y = ggplot2::element_blank(),
      #axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right"#, legend.text.align = 1
    ) +
    # ggplot2::annotate( # [NOTE] Damon removed this because it didn't seem to align with the diagonal.
    #   "raster", x=seq(nD), y=rev(seq(nD)), fill=divColor
    # ) +
    # four edges
    ggplot2::annotate("rect", xmin=.5-dw, xmax=.5,       ymin=.5-dw,            ymax=.5+nD+lim_expand+dw, fill=divColor) +
    ggplot2::annotate("rect", xmin=.5+nD, xmax=.5+nD+dw, ymin=.5-dw,            ymax=.5+nD+lim_expand+dw, fill=divColor) +
    ggplot2::annotate("rect", ymin=.5-dw, ymax=.5,       xmin=.5-lim_expand-dw, xmax=.5+nD+dw,            fill=divColor) +
    ggplot2::annotate("rect", ymin=.5+nD, ymax=.5+nD+dw, xmin=.5-lim_expand-dw, xmax=.5+nD+dw,            fill=divColor)

  if (use_groups) {
    # brain region label color bar
    plt <- plt + ggplot2::annotate(
      "rect", xmin=.5-lim_expand-dw, xmax=.5, ymin=.5+nD-seq(nD), ymax=.5+nD-seq(nD)+1,
      fill=group_cols[gdiv_cols]
    ) + ggplot2::annotate(
      "rect", ymin=.5+nD, ymax=.5+nD+lim_expand+dw, xmin=seq(nD)-.5, xmax=seq(nD)+.5,
      fill=group_cols[gdiv_cols]
    )
    xmin <- xmax <- NULL # to avoid missing var warnings in R package dev checks
    plt <- plt +
      # dividers
      ggplot2::geom_rect(data=gg_pdv, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=0, ymax=.5+nD+lim_expand+dw), fill=divColor) +
      ggplot2::geom_rect(data=gg_pdv, ggplot2::aes(ymin=nD+1-xmin, ymax=nD+1-xmax, xmin=.5-lim_expand-dw, xmax=.5+nD), fill=divColor) +
      # separate brain region label color bar
      ggplot2::annotate("rect", ymin=.5-dw, ymax=.5+nD+lim_expand+dw, xmin=.5-dw-.3, xmax=.5, fill="white") +
      ggplot2::annotate("rect", ymin=nD+.5, ymax=nD+.5+dw+.3, xmin=-.5-dw, xmax=.5+nD+dw, fill="white")
  }

  plt
}