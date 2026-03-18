# -------------------------------------------------------------------------
# ordinal_effect_measures
#
# Compute the causal estimands
#   eta   = P{Y(1) >  Y(0)}
#   tau   = P{Y(1) >= Y(0)}
# and their corresponding values under the working independence assumption
# with the same marginals:
#   eta_I = sum_j sum_{k<j} P(Y(1)=j) P(Y(0)=k)
#   tau_I = sum_j sum_{k<=j} P(Y(1)=j) P(Y(0)=k)
#
# INPUT
#   mat : numeric matrix representing a joint distribution or counts.
#
# BEHAVIOR
#   • If the entries do not sum to 1, the matrix is normalized.
#   • If the entries neither sum to 1 nor are all integers, the function
#     still normalizes the matrix but emits a warning.
#
# IMPORTANT: ORDER OF LEVELS
#   The matrix is assumed to be supplied in the natural order
#
#        rows: Y(1) = 1,2,...,J
#        cols: Y(0) = 1,2,...,J
#
#   This is the same interpretation used by the LaTeX table function.
#   Unlike the table function, this function does NOT reverse rows for
#   display, since it performs only numerical calculations.
#
# OUTPUT
#   A named list with elements:
#     eta, tau, eta_I, tau_I
#
# -------------------------------------------------------------------------

ordinal_effect_measures <- function(mat, normalize_tol = 1e-10) {
  if (!is.matrix(mat)) {
    stop("'mat' must be a matrix.")
  }
  
  if (!is.numeric(mat)) {
    stop("'mat' must be numeric.")
  }
  
  if (any(mat < 0, na.rm = TRUE)) {
    stop("'mat' must have nonnegative entries.")
  }
  
  s <- sum(mat)
  if (s <= 0) {
    stop("The matrix entries must sum to a positive value.")
  }
  
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
    stop("'mat' must be square, since Y(0) and Y(1) are assumed to have the same ordinal levels.")
  }
  
  J <- nr
  
  row_marg <- rowSums(mat)  # P(Y(1)=j)
  col_marg <- colSums(mat)  # P(Y(0)=k)
  
  eta <- 0
  tau <- 0
  
  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      if (j > k) {
        eta <- eta + mat[j, k]
      }
      if (j >= k) {
        tau <- tau + mat[j, k]
      }
    }
  }
  
  eta_I <- 0
  tau_I <- 0
  
  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      pij_ind <- row_marg[j] * col_marg[k]
      if (j > k) {
        eta_I <- eta_I + pij_ind
      }
      if (j >= k) {
        tau_I <- tau_I + pij_ind
      }
    }
  }
  
  list(
    eta = eta,
    tau = tau,
    eta_I = eta_I,
    tau_I = tau_I
  )
}