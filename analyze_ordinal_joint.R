# ------------------------------------------------------------
# Ordinal potential outcomes: joint-distribution analysis
#
# Input:
#   P : J x J joint probability matrix
#       rows    = Y(1)
#       columns = Y(0)
#
# The code calculates:
#
# 1. Local DTD conditions
#    - indices j where open-tail DTD holds:
#        P{Y(1) >  j | Y(0) = j} >= P{Y(1) >  j}
#    - indices j where closed-tail DTD holds:
#        P{Y(1) >= j | Y(0) = j} >= P{Y(1) >= j}
#
# 2. Estimands and bounds
#    - eta = P{Y(1) >  Y(0)}
#    - tau = P{Y(1) >= Y(0)}
#    - unrestricted lower and upper bounds
#    - independence quantities eta_I and tau_I
#
# 3. Local-DTD lower bounds
#    For each cutoff d = 1,...,J, the code calculates the
#    lower bounds obtained by assuming local DTD holds at
#    all indices k >= d.
#
#    cutoff = 1 therefore corresponds to global DTD and
#    gives the independence-based lower bounds.
#
#    cutoff = J means that local DTD is assumed only at J.
#
#    An additional row, cutoff = J + 1, corresponds to no
#    local-DTD assumption and therefore gives the unrestricted
#    lower bounds.
#
#    Open-tail local DTD is used to improve the lower bound
#    for eta.
#
#    Closed-tail local DTD is used to improve the lower bound
#    for tau.
#
#    The output also contains the individual candidate bounds
#    before maximizing over j.
#
# Matrix convention:
#   P[r, c] = P{Y(1) = r, Y(0) = c}
# ------------------------------------------------------------


analyze_ordinal_joint <- function(P, tol = 1e-12) {
  
  ## ------------------------------------------------------------
  ## Checks
  ## ------------------------------------------------------------
  
  P <- as.matrix(P)
  
  stopifnot(
    nrow(P) == ncol(P),
    all(P >= -tol),
    abs(sum(P) - 1) < 1e-8
  )
  
  P[P < 0 & P > -tol] <- 0
  
  J <- nrow(P)
  
  ## rows = Y(1), columns = Y(0)
  p1 <- rowSums(P)
  p0 <- colSums(P)
  
  
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
  ## 1. Indices at which local DTD holds
  ## ------------------------------------------------------------
  
  open_cond <- closed_cond <- rep(NA_real_, J)
  
  for (j in seq_len(J)) {
    
    if (p0[j] > tol) {
      
      # P{Y(1) > j | Y(0) = j}
      open_cond[j] <-
        if (j < J) sum(P[(j + 1):J, j]) / p0[j] else 0
      
      # P{Y(1) >= j | Y(0) = j}
      closed_cond[j] <-
        sum(P[j:J, j]) / p0[j]
    }
  }
  
  D_open <- which(open_cond + tol >= p1_gt)
  D_closed <- which(closed_cond + tol >= p1_ge)
  
  
  ## ------------------------------------------------------------
  ## 2. eta and tau
  ## ------------------------------------------------------------
  
  rr <- row(P)   # Y(1)
  cc <- col(P)   # Y(0)
  
  eta <- sum(P[rr > cc])
  tau <- sum(P[rr >= cc])
  
  
  ## ------------------------------------------------------------
  ## Unrestricted candidate bounds
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
  ## Independence quantities
  ## ------------------------------------------------------------
  
  P_ind <- outer(p1, p0)
  
  eta_I <- sum(P_ind[rr > cc])
  tau_I <- sum(P_ind[rr >= cc])
  
  
  ## ------------------------------------------------------------
  ## 3. Local-DTD lower bounds
  ##
  ## cutoff = d means assume local DTD at all k >= d.
  ##
  ## cutoff = J + 1 means no local-DTD assumption.
  ## ------------------------------------------------------------
  
  cutoffs <- 1:(J + 1)
  n_cut <- length(cutoffs)
  
  local_open_L <- rep(NA_real_, n_cut)
  local_closed_L <- rep(NA_real_, n_cut)
  
  open_candidate_matrix <- matrix(
    NA_real_,
    nrow = n_cut,
    ncol = J,
    dimnames = list(
      assumption_cutoff = cutoffs,
      bound_candidate = seq_len(J)
    )
  )
  
  closed_candidate_matrix <- matrix(
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
      
      open_candidate_matrix[ii, j] <-
        Delta[j] + extra
    }
    
    local_open_L[ii] <- max(open_candidate_matrix[ii, ])
    
    
    ## ----- closed-tail local DTD lower bound for tau -----
    
    for (j in seq_len(J)) {
      
      K <- D[D > j]
      
      extra <- if (length(K)) {
        sum(p1_ge[K] * p0[K])
      } else {
        0
      }
      
      closed_candidate_matrix[ii, j] <-
        p0[j] + Delta[j] + extra
    }
    
    local_closed_L[ii] <- max(closed_candidate_matrix[ii, ])
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
    
    DTD = list(
      open_indices = D_open,
      closed_indices = D_closed,
      open_conditional = open_cond,
      open_marginal = p1_gt,
      closed_conditional = closed_cond,
      closed_marginal = p1_ge
    ),
    
    estimands = c(
      eta = eta,
      tau = tau
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
    
    local_DTD_lower_bounds = data.frame(
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
      eta_open = local_open_L,
      tau_closed = local_closed_L
    ),
    
    local_DTD_candidates = list(
      eta_open = open_candidate_matrix,
      tau_closed = closed_candidate_matrix
    )
  )
}