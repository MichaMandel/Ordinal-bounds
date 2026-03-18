# -------------------------------------------------------------------------
# latex_joint_independence_table
#
# Generate LaTeX code for a table showing:
#   (1) a joint probability matrix,
#   (2) the corresponding independence distribution with the same marginals,
#   (3) row and column marginal probabilities.
#
# The layout matches the style used in ordinal causal inference tables,
# with two side-by-side blocks ("joint" and "independence") and margins.
#
# INPUT
#   mat : numeric matrix representing a joint distribution or counts.
#
# BEHAVIOR
#   • If the entries of the matrix do not sum to 1, they are normalized so
#     the matrix represents a probability distribution.
#
#   • If the entries neither sum to 1 nor are all integers (counts),
#     the function still normalizes the matrix but emits a warning.
#
#   • Row marginals P(Y(1)=y) and column marginals P(Y(0)=y) are computed.
#
#   • The independence distribution with the same marginals is computed as
#         P_ind(i,j) = P(Y(1)=i) P(Y(0)=j)
#
# IMPORTANT: ORDER OF LEVELS
#   The matrix is assumed to be supplied in the natural order
#
#        rows: Y(1) = 1,2,...,J
#        cols: Y(0) = 1,2,...,J
#
#   However, the LaTeX table is printed with the rows reversed so that the
#   highest ordinal level appears first:
#
#        Y(1) = J, J-1, ..., 1
#
#   This ordering matches the conventional presentation of contingency
#   tables in ordinal data analysis (largest category at the top).
#
# OUTPUT
#   A character string containing LaTeX code for a tabular environment
#   using the booktabs style. The code can be printed with:
#
#        cat(latex_joint_independence_table(mat))
#
#   and pasted directly into a LaTeX document.
#
# OPTIONAL ARGUMENTS
#   row_var  : label for the row variable (default "Y(1)")
#   col_var  : label for the column variable (default "Y(0)")
#   digits   : number of decimal places for probabilities
#   gap      : horizontal space between the two tables
#
# -------------------------------------------------------------------------

latex_joint_independence_table <- function(
    mat,
    row_var = "Y(1)",
    col_var = "Y(0)",
    row_marg_label = "P\\{Y(1)=y\\}",
    col_marg_label = "P\\{Y(0)=y\\}",
    digits = 4,
    gap = "1cm",
    normalize_tol = 1e-10
) {
  if (!is.matrix(mat)) stop("'mat' must be a matrix.")
  if (!is.numeric(mat)) stop("'mat' must be numeric.")
  if (any(mat < 0, na.rm = TRUE)) stop("'mat' must have nonnegative entries.")
  
  s <- sum(mat)
  if (s <= 0) stop("The matrix entries must sum to a positive value.")
  
  is_integerish <- function(x, tol = 1e-10) {
    all(abs(x - round(x)) < tol)
  }
  
  sums_to_one <- abs(s - 1) < normalize_tol
  
  if (!sums_to_one) {
    if (!is_integerish(mat)) {
      warning("Matrix entries do not sum to 1 and are not all integers; standardizing to sum to 1.")
    }
    mat <- mat / s
  }
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  if (nr != nc) {
    warning("Matrix is not square. The function will still create the table, but the layout is most natural for square matrices.")
  }
  
  row_marg <- rowSums(mat)
  col_marg <- colSums(mat)
  indep <- outer(row_marg, col_marg)
  
  fmt_num <- function(x, digits) {
    out <- formatC(x, format = "f", digits = digits)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    out
  }
  
  row_labels <- rev(seq_len(nr))
  col_labels <- seq_len(nc)
  
  mat_disp <- mat[rev(seq_len(nr)), , drop = FALSE]
  indep_disp <- indep[rev(seq_len(nr)), , drop = FALSE]
  row_marg_disp <- row_marg[rev(seq_len(nr))]
  
  colspec <- paste0(
    "c c ",
    paste(rep("c", nc), collapse = ""),
    " @{\\hspace{", gap, "}} ",
    paste(rep("c", nc), collapse = ""),
    " c"
  )
  
  lines <- character()
  lines <- c(lines, paste0("\\begin{tabular}{", colspec, "}"))
  lines <- c(lines, "\\toprule")
  lines <- c(
    lines,
    paste0(
      " &  & \\multicolumn{", nc, "}{c}{joint} & ",
      "\\multicolumn{", nc, "}{c}{independence} & \\\\"
    )
  )
  lines <- c(
    lines,
    paste0(
      "\\cmidrule(lr){3-", 2 + nc, "} ",
      "\\cmidrule(lr){", 3 + nc, "-", 2 + 2 * nc, "}"
    )
  )
  
  mid1 <- rep("", nc)
  mid1[ceiling(nc / 2)] <- paste0("{$", col_var, "$}")
  mid2 <- rep("", nc)
  mid2[ceiling(nc / 2)] <- paste0("{$", col_var, "$}")
  
  lines <- c(
    lines,
    paste0(
      " &  & ",
      paste(mid1, collapse = " & "),
      " & ",
      paste(mid2, collapse = " & "),
      " & \\\\"
    )
  )
  
  lines <- c(
    lines,
    paste0(
      " &  & ",
      paste(col_labels, collapse = " & "),
      " & ",
      paste(col_labels, collapse = " & "),
      " & $", row_marg_label, "$ \\\\"
    )
  )
  
  lines <- c(lines, "\\midrule")
  
  for (i in seq_len(nr)) {
    left_label <- if (i == ceiling(nr / 2)) paste0("$", row_var, "$") else ""
    row_joint <- fmt_num(mat_disp[i, ], digits)
    row_indep <- fmt_num(indep_disp[i, ], digits)
    row_marg_i <- fmt_num(row_marg_disp[i], digits)
    
    lines <- c(
      lines,
      paste0(
        left_label, " & ", row_labels[i], " & ",
        paste(row_joint, collapse = " & "), " & ",
        paste(row_indep, collapse = " & "), " & ",
        row_marg_i, " \\\\"
      )
    )
  }
  
  lines <- c(lines, "\\midrule")
  lines <- c(
    lines,
    paste0(
      " & $", col_marg_label, "$ & ",
      paste(fmt_num(col_marg, digits), collapse = " & "),
      " & ",
      paste(fmt_num(col_marg, digits), collapse = " & "),
      " & 1 \\\\"
    )
  )
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}")
  
  paste(lines, collapse = "\n")
}


