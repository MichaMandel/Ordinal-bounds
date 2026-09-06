# ------------------------------------------------------------
# Ordinal potential outcomes: bounds from marginal distributions
#
# Inputs:
#   p1 : marginal probabilities of Y(1), length J
#   p0 : marginal probabilities of Y(0), length J
#
# The code calculates quantities identified from the marginals
# or obtained under specified local-DTD assumptions.
#
# 1. Unrestricted bounds
#    - lower and upper bounds for
#        eta = P{Y(1) >  Y(0)}
#        tau = P{Y(1) >= Y(0)}
#    - all candidate lower and upper bounds before taking
#      the maximum/minimum over j
#
# 2. Independence quantities
#    - eta_I and tau_I, obtained under the working assumption
#      that Y(1) and Y(0) are independent
#
# 3. Local-DTD bounds
#    For each cutoff d = 1,...,J, the code calculates improved
#    lower and upper bounds obtained by assuming local DTD holds
#    at all indices k >= d.
#
#    cutoff = 1 therefore corresponds to global DTD.
#
#    cutoff = J means that local DTD is assumed only at J.
#
#    An additional row, cutoff = J + 1, corresponds to no
#    local-DTD assumption and therefore gives the unrestricted
#    bounds.
#
#    Open-tail local DTD is used to improve
#      - the lower bound for eta
#      - the upper bound for tau
#
#    Closed-tail local DTD is used to improve
#      - the lower bound for tau
#      - the upper bound for eta
#
#    The output also contains the individual candidate bounds
#    before maximizing/minimizing over j.
#
# Important:
#   These bounds depend only on the marginals p1 and p0.
#   No joint distribution of (Y(1), Y(0)) is required.
#
# Convention:
#   Categories are 1,...,J.
# ------------------------------------------------------------


analyze_ordinal_marginals <- function(p1, p0, tol = 1e-12) {
  
  ## ------------------------------------------------------------
  ## Checks
  ## ------------------------------------------------------------
  
  stopifnot(
    length(p1) == length(p0),
    all(p1 >= -tol),
    all(p0 >= -tol),
    abs(sum(p1) - 1) < 1e-8,
    abs(sum(p0) - 1) < 1e-8
  )
  
  p1[p1 < 0 & p1 > -tol] <- 0
  p0[p0 < 0 & p0 > -tol] <- 0
  
  J <- length(p1)
  
  
  ## ------------------------------------------------------------
  ## Marginal tail probabilities
  ## ------------------------------------------------------------
  
  # P{Y(1) >= j}, P{Y(0) >= j}
  p1_ge <- rev(cumsum(rev(p1)))
  p0_ge <- rev(cumsum(rev(p0)))
  
  # P{Y(1) > j}
  p1_gt <- c(p1_ge[-1], 0)
  
  # Delta_j = P{Y(1) >= j} - P{Y(0) >= j}
  Delta <- p1_ge - p0_ge
  
  
  ## ------------------------------------------------------------
  ## 1. Unrestricted candidate bounds
  ## ------------------------------------------------------------
  
  # eta
  eta_L_candidates <- Delta
  eta_U_candidates <- 1 + Delta - p1
  
  eta_L <- max(eta_L_candidates)
  eta_U <- min(eta_U_candidates)
  
  # tau
  tau_L_candidates <- p0 + Delta
  tau_U_candidates <- 1 + Delta
  
  tau_L <- max(tau_L_candidates)
  tau_U <- min(tau_U_candidates)
  
  
  ## ------------------------------------------------------------
  ## 2. Independence quantities
  ## ------------------------------------------------------------
  
  # P_ind[r, c] = P{Y(1) = r} P{Y(0) = c}
  P_ind <- outer(p1, p0)
  
  rr <- row(P_ind)
  cc <- col(P_ind)
  
  eta_I <- sum(P_ind[rr > cc])
  tau_I <- sum(P_ind[rr >= cc])
  
  
  ## ------------------------------------------------------------
  ## 3. Local-DTD bounds
  ##
  ## cutoff = d means assume local DTD at all k >= d.
  ##
  ## cutoff = J + 1 means no local-DTD assumption.
  ## ------------------------------------------------------------
  
  cutoffs <- 1:(J + 1)
  n_cut <- length(cutoffs)
  
  local_open_L   <- rep(NA_real_, n_cut)
  local_closed_L <- rep(NA_real_, n_cut)
  local_open_U   <- rep(NA_real_, n_cut)
  local_closed_U <- rep(NA_real_, n_cut)
  
  open_L_candidate_matrix <- matrix(
    NA_real_,
    nrow = n_cut,
    ncol = J,
    dimnames = list(
      assumption_cutoff = cutoffs,
      bound_candidate = seq_len(J)
    )
  )
  
  closed_L_candidate_matrix <- matrix(
    NA_real_,
    nrow = n_cut,
    ncol = J,
    dimnames = list(
      assumption_cutoff = cutoffs,
      bound_candidate = seq_len(J)
    )
  )
  
  open_U_candidate_matrix <- matrix(
    NA_real_,
    nrow = n_cut,
    ncol = J,
    dimnames = list(
      assumption_cutoff = cutoffs,
      bound_candidate = seq_len(J)
    )
  )
  
  closed_U_candidate_matrix <- matrix(
    NA_real_,
    nrow = n_cut,
    ncol = J,
    dimnames = list(
      assumption_cutoff = cutoffs,
      bound_candidate = seq_len(J)
    )
  )
  
  
  for (ii in seq_along(cutoffs)) {
    
    d <- cutoffs[ii]
    
    # Assume local DTD at all k >= d
    D <- if (d <= J) d:J else integer(0)
    
    
    ## ----- open-tail local DTD lower bound for eta -----
    
    for (j in seq_len(J)) {
      
      K <- D[D >= j]
      
      extra <- if (length(K)) {
        sum(p1_gt[K] * p0[K])
      } else {
        0
      }
      
      open_L_candidate_matrix[ii, j] <-
        Delta[j] + extra
    }
    
    local_open_L[ii] <- max(open_L_candidate_matrix[ii, ])
    
    
    ## ----- closed-tail local DTD lower bound for tau -----
    
    for (j in seq_len(J)) {
      
      K <- D[D > j]
      
      extra <- if (length(K)) {
        sum(p1_ge[K] * p0[K])
      } else {
        0
      }
      
      closed_L_candidate_matrix[ii, j] <-
        p0[j] + Delta[j] + extra
    }
    
    local_closed_L[ii] <- max(closed_L_candidate_matrix[ii, ])
    
    
    ## ----- open-tail local DTD upper bound for tau -----
    
    for (j in seq_len(J)) {
      
      improvement <- if (j > 1 && (j - 1) %in% D) {
        p1_gt[j - 1] * p0[j - 1]
      } else {
        0
      }
      
      open_U_candidate_matrix[ii, j] <-
        1 + Delta[j] - improvement
    }
    
    local_open_U[ii] <- min(open_U_candidate_matrix[ii, ])
    
    
    ## ----- closed-tail local DTD upper bound for eta -----
    
    for (j in seq_len(J)) {
      
      improvement <- if (j < J && J %in% D) {
        p1[J] * p0[J]
      } else {
        0
      }
      
      closed_U_candidate_matrix[ii, j] <-
        1 + Delta[j] - p1[j] - improvement
    }
    
    local_closed_U[ii] <- min(closed_U_candidate_matrix[ii, ])
  }
  
  
  ## ------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------
  
  list(
    
    marginals = list(
      p1 = p1,
      p0 = p0,
      p1_ge = p1_ge,
      p1_gt = p1_gt,
      p0_ge = p0_ge,
      Delta = Delta
    ),
    
    unrestricted_bounds = list(
      eta = c(lower = eta_L, upper = eta_U),
      tau = c(lower = tau_L, upper = tau_U),
      eta_lower_candidates = eta_L_candidates,
      eta_upper_candidates = eta_U_candidates,
      tau_lower_candidates = tau_L_candidates,
      tau_upper_candidates = tau_U_candidates
    ),
    
    independence = c(
      eta_I = eta_I,
      tau_I = tau_I
    ),
    
    local_DTD_bounds = data.frame(
      cutoff = cutoffs,
      assume_DTD_for = sapply(
        cutoffs,
        function(d) {
          if (d <= J) {
            paste(d:J, collapse = ",")
          } else {
            "none"
          }
        }
      ),
      eta_lower_open = local_open_L,
      eta_upper_closed = local_closed_U,
      tau_lower_closed = local_closed_L,
      tau_upper_open = local_open_U
    ),
    
    local_DTD_candidates = list(
      eta_lower_open = open_L_candidate_matrix,
      eta_upper_closed = closed_U_candidate_matrix,
      tau_lower_closed = closed_L_candidate_matrix,
      tau_upper_open = open_U_candidate_matrix
    )
  )
}
