# Colored probability matrices for ordinal causal effects
#
# Columns represent Y(0), rows represent Y(1), and cell (x, y) represents
# p[x,y] = P{Y(0) = x, Y(1) = y}.
#
# Colors follow the paper's convention:
#   green  = add the cell probability
#   orange = subtract the cell probability
#   white  = coefficient zero
# A "x2" mark means that the cell probability is added twice.

col_plus  <- "#009E73"  # color-blind-accessible green
col_minus <- "#D55E00"  # color-blind-accessible orange


# Return the coefficient attached to every cell of the J by J matrix.
# This helper is useful if the cell definitions need to be inspected directly.
colored_matrix_coefficients <- function(
    J, j, quantity, term, difference = FALSE
) {
  if (length(J) != 1L || !is.numeric(J) || is.na(J) ||
      J < 2 || J != as.integer(J)) {
    stop("J must be a single integer greater than or equal to 2.")
  }
  J <- as.integer(J)

  if (length(j) != 1L || !is.numeric(j) || is.na(j) ||
      j < 1 || j > J || j != as.integer(j)) {
    stop("j must be a single integer between 1 and J.")
  }
  j <- as.integer(j)

  if (length(quantity) != 1L || !is.character(quantity)) {
    stop("quantity must be one of 'delta', 'tau', or 'eta'.")
  }
  quantity <- tolower(quantity)
  if (!quantity %in% c("delta", "tau", "eta")) {
    stop("quantity must be one of 'delta', 'tau', or 'eta'.")
  }

  if (length(term) != 1L || !is.character(term)) {
    stop("term must be one of 'L', 'U', or 'E'.")
  }
  term <- toupper(term)
  if (!term %in% c("L", "U", "E")) {
    stop("term must be one of 'L', 'U', or 'E'.")
  }
  if (quantity == "delta" && term != "E") {
    stop("delta_j is exact; use term = 'E'.")
  }
  if (length(difference) != 1L || !is.logical(difference) ||
      is.na(difference)) {
    stop("difference must be either TRUE or FALSE.")
  }
  if (difference && term == "E") {
    stop("To plot a bound minus the exact value, term must be 'L' or 'U'.")
  }
  if (difference && quantity == "delta") {
    stop("delta_j has no lower or upper bound in this construction.")
  }

  cells <- expand.grid(x = seq_len(J), y = seq_len(J))
  x <- cells$x
  y <- cells$y

  # Delta_j = P{Y(1) >= j} - P{Y(0) >= j}.
  delta_coef <- as.integer(y >= j) - as.integer(x >= j)

  coef <- switch(
    paste(quantity, term, sep = "_"),

    # Exact quantities
    delta_E = delta_coef,
    tau_E   = as.integer(y >= x),
    eta_E   = as.integer(y > x),

    # Candidate lower and upper terms indexed by j:
    # tau_L(j) = P{Y(0) = j} + Delta_j
    tau_L = as.integer(x == j) + delta_coef,
    # tau_U(j) = 1 + Delta_j
    tau_U = 1L + delta_coef,
    # eta_L(j) = Delta_j
    eta_L = delta_coef,
    # eta_U(j) = 1 + Delta_j - P{Y(1) = j}
    eta_U = 1L + delta_coef - as.integer(y == j),

    stop("This combination of quantity and term is not defined.")
  )

  if (difference) {
    exact_coef <- switch(
      quantity,
      tau = as.integer(y >= x),
      eta = as.integer(y > x)
    )
    coef <- coef - exact_coef
  }

  data.frame(cells, coefficient = coef)
}


# Draw one colored matrix.
#
# Required arguments:
#   J        number of ordinal categories
#   j        cutoff index (1, ..., J)
#   quantity "delta", "tau", or "eta"
#   term     "E" (exact), "L" (candidate lower term), or
#            "U" (candidate upper term)
#   difference  if TRUE, plot the candidate bound minus the exact quantity
#
# Optional arguments:
#   info                if TRUE, add an automatically generated information title
#   show_probabilities  add p[x,y] labels to the cells
#   main                custom plot title; overrides the automatic title
#   cex                 scaling for axes and cell annotations
#
# Supported combinations are:
#   delta/E; tau/E, tau/L, tau/U; eta/E, eta/L, eta/U.
# For L and U, the plot is the candidate term corresponding to the supplied j;
# the overall sharp bound is obtained by maximizing/minimizing across j.
plot_colored_ordinal_matrix <- function(
    J, j, quantity, term,
    difference = FALSE,
    info = TRUE,
    show_probabilities = FALSE,
    main = NULL,
    cex = 1
) {
  cells <- colored_matrix_coefficients(
    J, j, quantity, term, difference = difference
  )
  J <- as.integer(J)
  j <- as.integer(j)
  quantity <- tolower(quantity)
  term <- toupper(term)

  if (length(info) != 1L || !is.logical(info) || is.na(info)) {
    stop("info must be either TRUE or FALSE.")
  }

  # Construct a two-line title unless a custom title was supplied. The first
  # line identifies j and the parameter; the second identifies what is drawn.
  if (info && is.null(main)) {
    parameter_symbol <- switch(
      quantity,
      delta = bquote(Delta[.(j)]),
      tau = quote(tau),
      eta = quote(eta)
    )

    term_description <- switch(
      term,
      E = "exact value",
      L = "candidate lower bound",
      U = "candidate upper bound"
    )
    if (difference) {
      term_description <- paste("difference:", term_description)
    }

    main <- bquote(atop(
      italic(j) == .(j) * "," ~~ "parameter:" ~~ .(parameter_symbol),
      .(term_description)
    ))
  }

  # A square plotting region keeps the axes adjacent to the square matrix even
  # when the graphics device itself is rectangular.
  old_par <- par(xaxs = "i", yaxs = "i", pty = "s")
  on.exit(par(old_par), add = TRUE)

  plot.new()
  plot.window(xlim = c(0.5, J + 0.5), ylim = c(0.5, J + 0.5))

  fill <- rep("white", nrow(cells))
  fill[cells$coefficient > 0] <- col_plus
  fill[cells$coefficient < 0] <- col_minus

  rect(
    cells$x - 0.5, cells$y - 0.5,
    cells$x + 0.5, cells$y + 0.5,
    col = fill, border = "#808080", lwd = 0.8
  )
  rect(0.5, 0.5, J + 0.5, J + 0.5, border = "black", lwd = 1.2)

  axis(1, at = seq_len(J), labels = seq_len(J), tick = FALSE,
       line = 0, cex.axis = cex)
  axis(2, at = seq_len(J), labels = seq_len(J), tick = FALSE,
       line = 0, las = 1, cex.axis = cex)
  title(main = main, cex.main = cex)
  mtext(expression(italic(Y)^{(0)}), side = 1, line = 1.8,
        cex = cex)
  mtext(expression(italic(Y)^{(1)}), side = 2, line = 1.8,
        cex = cex)

  if (show_probabilities) {
    labs <- as.expression(lapply(seq_len(nrow(cells)), function(k) {
      bquote(italic(p)[.(cells$x[k]) * "," * .(cells$y[k])])
    }))
    text(cells$x, cells$y, labels = labs, cex = 0.72 * cex)
  }

  doubled <- cells$coefficient == 2L
  if (any(doubled)) {
    text(cells$x[doubled], cells$y[doubled], labels = "×2",
         cex = 0.9 * cex)
  }

  invisible(cells)
}


# Examples ---------------------------------------------------------------
#
# Exact Delta_3:
# plot_colored_ordinal_matrix(6, 3, "delta", "E",
#                             show_probabilities = TRUE)
#
# Candidate lower and upper terms for tau at j = 3:
# plot_colored_ordinal_matrix(6, 3, "tau", "L")
# plot_colored_ordinal_matrix(6, 3, "tau", "U")
#
# Exact eta and its candidate upper term at j = 3:
# plot_colored_ordinal_matrix(6, 3, "eta", "E")
# plot_colored_ordinal_matrix(6, 3, "eta", "U")
#
# Difference between a candidate bound and the exact quantity:
# plot_colored_ordinal_matrix(6, 3, "tau", "L", difference = TRUE)
# plot_colored_ordinal_matrix(6, 3, "tau", "U", difference = TRUE)
# plot_colored_ordinal_matrix(6, 3, "eta", "L", difference = TRUE)
# plot_colored_ordinal_matrix(6, 3, "eta", "U", difference = TRUE)
#
# Suppress the automatic information title:
# plot_colored_ordinal_matrix(6, 3, "tau", "U", info = FALSE)
#
# Save a plot to PDF, for example:
# pdf("tau_upper_j3.pdf", width = 4, height = 4)
# plot_colored_ordinal_matrix(6, 3, "tau", "U")
# dev.off()
