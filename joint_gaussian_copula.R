# ------------------------------------------------------------
# Gaussian copula model for ordinal potential outcomes
#
# Purpose:
#   Construct a J x J joint probability matrix for
#   (Y(1), Y(0)) using a latent bivariate normal model
#   with correlation parameter rho.
#
# Inputs:
#   p1  : marginal probabilities of Y(1), length J
#   p0  : marginal probabilities of Y(0), length J
#   rho : latent Gaussian correlation, with -1 < rho < 1
#
# Method:
#   - Convert the marginal cumulative probabilities into
#     standard-normal cutpoints.
#   - Assume latent variables (Z1, Z0) follow a bivariate
#     standard normal distribution with correlation rho.
#   - Integrate the bivariate normal density over the
#     corresponding rectangles to obtain the joint cell
#     probabilities.
#
# Output:
#   P : J x J joint probability matrix with
#       rows    = Y(1)
#       columns = Y(0)
#
# Thus:
#   P[j, k] = P{Y(1) = j, Y(0) = k}
#
# For fixed marginals p1 and p0, rho controls the dependence
# structure while leaving both marginal distributions unchanged.
#
# Requires:
#   package mvtnorm
# ------------------------------------------------------------

joint_gaussian_copula <- function(p1, p0, rho) {
  stopifnot(
    length(p1) == length(p0),
    all(p1 >= 0),
    all(p0 >= 0),
    abs(sum(p1) - 1) < 1e-8,
    abs(sum(p0) - 1) < 1e-8,
    abs(rho) < 1
  )
  
  J <- length(p1)
  
  # Latent cutpoints, including -Inf and Inf
  a1 <- qnorm(c(0, cumsum(p1)))
  a0 <- qnorm(c(0, cumsum(p0)))
  
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  P <- matrix(0, nrow = J, ncol = J)
  
  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      P[j, k] <- mvtnorm::pmvnorm(
        lower = c(a1[j],     a0[k]),
        upper = c(a1[j + 1], a0[k + 1]),
        mean  = c(0, 0),
        sigma = Sigma
      )[1]
    }
  }
  
  dimnames(P) <- list(
    `Y(1)` = seq_len(J),
    `Y(0)` = seq_len(J)
  )
  
  P
}