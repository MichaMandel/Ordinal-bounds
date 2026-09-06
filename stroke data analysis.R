# analysis of stroke data from NEJM 2015
# data of n=500 patients by intention-to-treat
# data taken from Figure 1 of the paper
# outcome is the modified rankin scale score after 90 days
# score is from 0 to 6:
# 0 - no symptoms, 
# 1 - no clinically significant disability, 
# 2 - slight disability (patient is able to look after own affairs 
#     without assistance but is unable to carry 
#     out all previous activities), 
# 3 - moderate disability (patient requires some 
#     help but is able to walk unassisted),
# 4 - moderately severe disability (patient is 
#     unable to attend to bodily needs without 
#     assistance and unable to walk unassisted), 
# 5 - severe disability (patient requires constant
#     nursing care and attention), 
# 6 - death.
#
# treatment - standard treatment + Intraarterial 
#     treatment consisted of arterial catheterization
#     with a microcatheter to the level of
#     occlusion and delivery of a thrombolytic agent,
#     mechanical thrombectomy, or both.
# control - standard treatment (medication, maybe using IV)

# Data from Table S5 of the supplementary appendix
# combine 0 and 1  because of small n's
score <- 1:6
trt <- c(27,49,43,52,13,49)
ctrl <- c(16,35,44,81,32,59)

# reverse the order so larger index = better outcome
trt.ord <- trt[6:1]
ctrl.ord <- ctrl[6:1]

# call the function for the analysis
source("analyze_ordinal_marginals.R")

p.trt <- trt.ord/sum(trt.ord)
p.ctrl <- ctrl.ord/sum(ctrl.ord)

res <- analyze_ordinal_marginals(p1 = p.trt , 
                                 p0 = p.ctrl)

round(p.trt,3)
round(p.ctrl,3)
round(res$marginals$p1_ge,3)
round(res$marginals$p0_ge,3)
round(res$marginals$Delta,3)

round(res$unrestricted_bounds$eta,3)
round(res$unrestricted_bounds$tau,3)
round(res$independence,3)

res$local_DTD_bounds

# table for the latex file 

tab <- res$local_DTD_bounds[
  ,
  c(
    "assume_DTD_for",
    "eta_lower_open",
    "eta_upper_closed",
    "tau_lower_closed",
    "tau_upper_open"
  )
]

latex_tab <- knitr::kable(
  tab,
  format = "latex",
  digits = 3,
  booktabs = TRUE,
  escape = FALSE,
  col.names = c(
    "DTD set$^{*}$",
    "$\\tilde{\\eta}_L$",
    "$\\tilde{\\eta}_U$",
    "$\\tilde{\\tau}_L$",
    "$\\tilde{\\tau}_U$"
  )
)

latex_tab <- gsub("\\\\addlinespace\\n?", "", latex_tab)

cat(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  latex_tab,
  "\n\\begin{minipage}{0.95\\textwidth}\n",
  "\\footnotesize\n",
  "$^{*}$ The indicated set is assumed to satisfy closed-tail DTD for ",
  "$\\tilde{\\tau}_L$ and $\\tilde{\\eta}_U$, and open-tail DTD for ",
  "$\\tilde{\\eta}_L$ and $\\tilde{\\tau}_U$.\n",
  "\\end{minipage}\n",
  "\\caption{Bounds under local DTD assumptions.}\n",
  "\\label{tab:local_DTD}\n",
  "\\end{table}\n",
  sep = ""
)
######################################################
#
#           SENSITIVITY ANALYSIS
#
####################################################

# call the function for the sensitivity analysis
source("analyze_ordinal_joint.R")
source("joint_gaussian_copula.R")

rho = 0.3
P.sens <- joint_gaussian_copula(p1 = p.trt,
                                p0 = p.ctrl,
                                rho = rho)
res.joint <- analyze_ordinal_joint(P.sens)

# check that we get the same results for bounds

round(res.joint$marginals$p1,3)
round(res.joint$marginals$p0,3)
round(res.joint$marginals$Delta,3)

round(res.joint$unrestricted_bounds$eta,3)
round(res.joint$unrestricted_bounds$tau,3)
round(res.joint$independence,3)

res.joint$local_DTD_lower_bounds

# run the analysis for a sequence of rho values
rhos <- round(seq(-0.9, 0.9, by = 0.1), 1)

results <- lapply(rhos, function(rho) {
  P <- joint_gaussian_copula(p.trt, p.ctrl, rho)
  res <- analyze_ordinal_joint(P)
  res$rho <- rho
  res
})

names(results) <- sprintf("rho_%+.1f", rhos)


# plotting eta and tau vs the lower bounds

library(ggplot2)

# Extract estimands across rho
plot_df <- data.frame(
  rho = rhos,
  eta = sapply(results, function(x) x$estimands["eta"]),
  tau = sapply(results, function(x) x$estimands["tau"])
)

# Bounds depend only on the marginals, so take them from one result
eta_bounds <- results[[1]]$local_DTD_lower_bounds$eta_open
tau_bounds <- results[[1]]$local_DTD_lower_bounds$tau_closed


p_eta <- ggplot(plot_df, aes(x = rho, y = eta)) +
  geom_hline(yintercept = eta_bounds, linetype = "dashed",
             linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(-0.9, 0.9, by = 0.3),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    x = expression(rho),
    y = expression(eta)
  ) +
  theme_classic(base_size = 12)

p_eta


p_tau <- ggplot(plot_df, aes(x = rho, y = tau)) +
  geom_hline(yintercept = tau_bounds, linetype = "dashed",
             linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(-0.9, 0.9, by = 0.3),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    x = expression(rho),
    y = expression(tau)
  ) +
  theme_classic(base_size = 12)

p_tau


ggsave(
  filename = file.path("figures", "eta_gaussian_copula.pdf"),
  plot = p_eta,
  width = 6,
  height = 4.5
)

ggsave(
  filename = file.path("figures", "tau_gaussian_copula.pdf"),
  plot = p_tau,
  width = 6,
  height = 4.5
)

#####################################
#####################################
# indices for which local DTD holds

DTD_table <- data.frame(
  rho = rhos,
  open_DTD = sapply(
    results,
    function(x) paste(x$DTD$open_indices, collapse = ",")
  ),
  closed_DTD = sapply(
    results,
    function(x) paste(x$DTD$closed_indices, collapse = ",")
  )
)

DTD_table
