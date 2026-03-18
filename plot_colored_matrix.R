#' Plot a K × K colored matrix grid
#'
#' Draws a square K-by-K grid with axis labels 1..K centered in each cell.
#' Selected cells can be filled with specified colors, and optional labels
#' (e.g., ×2) can be added on top of selected cells.
#'
#' @param K Integer. Number of rows and columns in the grid.
#' @param sets A list of lists. Each inner list must contain:
#'   \itemize{
#'     \item x: vector of column indices (1..K)
#'     \item y: vector of row indices (1..K)
#'     \item col: a color name (or vector of colors)
#'   }
#'   Each set defines a group of cells to be filled with the specified color.
#'
#' @param marks Optional list specifying cells where a label should be added.
#'   Must contain:
#'   \itemize{
#'     \item x: vector of column indices (1..K)
#'     \item y: vector of row indices (1..K)
#'   }
#'
#' @param mark_label Plotmath expression used as the label drawn in marked
#'   cells. Default is \eqn{\times 2}.
#'
#' @param label_cex Numeric scaling factor controlling the size of the mark
#'   labels. The default scales automatically with \code{K} so that labels
#'   remain inside cells when the grid becomes large.
#'
#' @param main Optional plot title.
#' @param xlab Label for the x-axis (can use \code{expression()} for math).
#' @param ylab Label for the y-axis (can use \code{expression()} for math).
#' @param grid_col Color of grid lines.
#' @param border_col Border color of filled cells (default \code{NA}).
#'
#' @details
#' Cells are indexed so that (x, y) corresponds to column x and row y,
#' with (1,1) at the bottom-left corner of the grid.
#'
#' If a cell appears in multiple sets, later sets overwrite earlier ones.
#' Labels specified via \code{marks} are drawn on top of the colored cells.
#'
#' @examples
#' x1 <- c(1,2,3)
#' y1 <- c(1,1,2)
#'
#' plot_colored_matrix(
#'   K = 7,
#'   sets = list(list(x = x1, y = y1, col = "red")),
#'   marks = list(x = c(2,3), y = c(1,2)),
#'   xlab = expression(Y^{(0)}),
#'   ylab = expression(Y^{(1)})
#' )


plot_colored_matrix <- function(
    K = 7,
    sets = list(),              # list of list(x=..., y=..., col=...)
    marks = NULL,               # adding labels to the cells
    mark_label = expression("\u00D7" * 2), #default label X2
    label_cex = min(1.3, 10/K),
    main = NULL, 
    xlab = NULL, 
    ylab = NULL,
    grid_col = "grey50",
    border_col = NA             # set to "grey40" if you want borders on filled cells
) {
  # Basic validation
  if (!is.list(sets)) stop("'sets' must be a list of sets.")
  
  op <- par(mar = c(4, 4, 2, 2) + 0.1, pty = "s", mgp = c(1.8, 0.4, 0))
  on.exit(par(op), add = TRUE)
  
  plot(NA,
       xlim = c(0, K), ylim = c(0, K),
       xaxs = "i", yaxs = "i",
       asp = 1, axes = FALSE,
       main = main, xlab = xlab, ylab = ylab)
  
  # Fill colored cells first (so grid lines draw on top)
  for (s in sets) {
    if (is.null(s$x) || is.null(s$y) || is.null(s$col)) {
      stop("Each set must be a list with elements x, y, col.")
    }
    x <- s$x; y <- s$y; col <- s$col
    if (length(x) != length(y)) stop("In each set, x and y must have same length.")
    
    # Filter to valid cells
    ok <- (x >= 1 & x <= K & y >= 1 & y <= K)
    x <- x[ok]; y <- y[ok]
    if (length(x) == 0) next
    
    # Draw rectangles for each (x,y): cell spans [x-1,x] × [y-1,y]
    rect(xleft = x - 1, ybottom = y - 1,
         xright = x,    ytop    = y,
         col = col, border = border_col)
  }
  
  # Grid on top
  abline(v = 0:K, h = 0:K, col = grid_col)
  
  # add labels (defualt ×2) markers
  if (!is.null(marks)) {
    
    n <- length(marks$x)
    
    text(
      marks$x - 0.5,
      marks$y - 0.5,
      labels = mark_label,
      font = 2,
      cex = label_cex
    )
    
  }
  
  # Axis labels at cell centers
  axis(1, at = (1:K) - 0.5, labels = 1:K, tick = FALSE)
  axis(2, at = (1:K) - 0.5, labels = 1:K, tick = FALSE)
  box()
}


