# tables for the paper and estimands


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("latex_joint_independence_table.R")
source("ordinal_effect_measures.R")

# PRD holds but the bounds do not
mat <- matrix(
  c(0.02, 0.04, 0.06,
    0.09, 0.07, 0.48,
    0.03, 0.02, 0.19),
  nrow = 3, byrow = TRUE
)

cat(latex_joint_independence_table(mat, digits = 4))

ordinal_effect_measures(mat)


# closed- but not open-tail holds, eta<eta_I

mat <- matrix(
  c(0.20, 0.06, 0.04,
    0.08, 0.22, 0.05,
    0.07, 0.18, 0.10),
  nrow = 3, byrow = TRUE
)

cat(latex_joint_independence_table(mat, digits = 4))
ordinal_effect_measures(mat)


# open- but not closed-tail holds, tau<tau_I

mat <- matrix(
c(0.05, 0.10, 0.10,
0.30, 0.05, 0.05,
0.20, 0.10, 0.05),
nrow = 3, byrow = TRUE
)

cat(latex_joint_independence_table(mat, digits = 4))

ordinal_effect_measures(mat)


# both conditions hold

mat <- matrix(
  c(0.16, 0.11, 0.12,
0.11, 0.10, 0.05,
0.13, 0.13, 0.09),
nrow = 3, byrow = TRUE
)

cat(latex_joint_independence_table(mat, digits = 4))
ordinal_effect_measures(mat)
