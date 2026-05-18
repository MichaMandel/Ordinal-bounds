# Improved lower bounds for eta and tau
#
# Usage:
#   improved_bounds(y0, y1, D_ot = integer(0), D_ct = integer(0))
#
# Arguments:
#   y0    A vector of probabilities or counts for the marginal distribution of Y(0).
#   y1    A vector of probabilities or counts for the marginal distribution of Y(1).
#         The vectors y0 and y1 must have the same length. The first entry is
#         assumed to correspond to outcome level 1, the second to level 2, and so on.
#         Higher levels correspond to better outcomes.
#
#   D_ot  A vector of indices at which the local open-tail DTD condition is assumed
#         to hold. These indices must be between 1 and J, where J = length(y0).
#
#   D_ct  A vector of indices at which the local closed-tail DTD condition is assumed
#         to hold. These indices must be between 1 and J.
#
# Defaults:
#   If D_ot and D_ct are not specified, they are taken to be empty. In this case,
#   the function returns the nonparametric sharp lower bounds.
#
#   If both D_ot and D_ct are specified as 1:J, the function returns the
#   independence-based lower bounds.
#
# Notes:
#   Inputs may be either probabilities or counts. Count vectors are automatically
#   normalized to probabilities. If an input vector does not sum to 1 and does not
#   appear to be a vector of counts, the function returns a warning.
#
# Output:
#   A list containing:
#     eta_lower  The lower bound for eta.
#     tau_lower  The lower bound for tau.
#     eta_terms  The candidate lower bounds for eta inside the maximum for each j.
#     tau_terms  The candidate lower bounds for tau inside the maximum for each j.
#     Delta      The vector Delta_j = P{Y(1) >= j} - P{Y(0) >= j}.
#     p0, p1     The normalized marginal distributions of Y(0) and Y(1).

improved_bounds <- function(y0, y1, D_ot = integer(0), D_ct = integer(0), tol = 1e-8) {
  
  # Basic checks
  if (length(y0) != length(y1)) {
    stop("y0 and y1 must have the same length.")
  }
  
  if (any(y0 < 0) || any(y1 < 0)) {
    stop("All entries of y0 and y1 must be nonnegative.")
  }
  
  J <- length(y0)
  
  if (any(D_ot < 1 | D_ot > J)) {
    stop("All indices in D_ot must be between 1 and length(y0).")
  }
  
  if (any(D_ct < 1 | D_ct > J)) {
    stop("All indices in D_ct must be between 1 and length(y0).")
  }
  
  # Check whether inputs look like probabilities or counts
  sum_y0 <- sum(y0)
  sum_y1 <- sum(y1)
  
  y0_is_prob <- abs(sum_y0 - 1) < tol
  y1_is_prob <- abs(sum_y1 - 1) < tol
  
  y0_is_count <- all(abs(y0 - round(y0)) < tol)
  y1_is_count <- all(abs(y1 - round(y1)) < tol)
  
  if (!y0_is_prob && !y0_is_count) {
    warning("y0 does not sum to 1 and does not appear to be a vector of counts.")
  }
  
  if (!y1_is_prob && !y1_is_count) {
    warning("y1 does not sum to 1 and does not appear to be a vector of counts.")
  }
  
  # Convert counts to probabilities if needed
  p0 <- y0 / sum_y0
  p1 <- y1 / sum_y1
  
  # Delta_j = P{Y(1) >= j} - P{Y(0) >= j}
  tail0_closed <- rev(cumsum(rev(p0)))  # P{Y(0) >= j}
  tail1_closed <- rev(cumsum(rev(p1)))  # P{Y(1) >= j}
  Delta <- tail1_closed - tail0_closed
  
  # P{Y(1) > k}
  tail1_open <- c(tail1_closed[-1], 0)
  
  eta_terms <- numeric(J)
  tau_terms <- numeric(J)
  
  for (j in seq_len(J)) {
    
    # eta lower bound term:
    # Delta_j + sum_{k >= j, k in D_ot} P{Y(1) > k} P{Y(0)=k}
    k_eta <- D_ot[D_ot >= j]
    
    eta_terms[j] <- Delta[j] +
      sum(tail1_open[k_eta] * p0[k_eta])
    
    # tau lower bound term:
    # P{Y(0)=j} + Delta_j + sum_{k > j, k in D_ct} P{Y(1) >= k} P{Y(0)=k}
    k_tau <- D_ct[D_ct > j]
    
    tau_terms[j] <- p0[j] + Delta[j] +
      sum(tail1_closed[k_tau] * p0[k_tau])
  }
  
  list(
    eta_lower = max(eta_terms),
    tau_lower = max(tau_terms),
    eta_terms = eta_terms,
    tau_terms = tau_terms,
    Delta = Delta,
    p0 = p0,
    p1 = p1
  )
}