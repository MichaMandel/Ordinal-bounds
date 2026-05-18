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

# Data from Table S5 of the supplemantery appendix
# combine 0 and 1  because of small n's
score <- 1:6
trt <- c(27,49,43,52,13,49)
ctrl <- c(16,35,44,81,32,59)
P.trt <- cumsum(trt)/sum(trt)
P.ctrl <- cumsum(ctrl)/sum(ctrl)
marginal.trt <- trt/sum(trt)
marginal.ctrl <- ctrl/sum(ctrl)
Delta <- P.trt-P.ctrl
round(rbind(score,trt,ctrl,marginal.trt,marginal.ctrl,
            P.trt,P.ctrl,Delta),3)
tau_L = max(marginal.ctrl+Delta)
tau_U = 1 + min(Delta)
eta_L = max(Delta)
eta_U = 1+min(Delta-marginal.trt)
round(cbind(tau_L,tau_U,eta_L,eta_U),3)


# bounds under independence 
joint = (outer(marginal.trt,marginal.ctrl,"*"))
tau_I = 0
eta_I = 0
for (i in 1:6) {
  for (j in i:6) {
    tau_I = tau_I + joint[i,j]
    if (j>i) {eta_I = eta_I + joint[i,j]}
  }
}
tau_I
eta_I

# the weighted effect with w the average of the probabilities
w = (marginal.ctrl[1:5]+marginal.trt[1:5])/2
w <- w/sum(w)
sum(w[1:5]*Delta[1:5])

# Using the function for lower bounds
source("improved lower bounds")

# (levels in reversed order: from worse to best)

#Model-free bounds 
improved_bounds(y0=ctrl[6:1], y1=trt[6:1], D_ot = integer(0), D_ct = integer(0), tol = 1e-8)

# assuming local DTD at 5 and 6
improved_bounds(y0=ctrl[6:1], y1=trt[6:1], D_ot = c(5,6), D_ct = c(5,6), tol = 1e-8)

# assuming local DTD at 4, 5 and 6
improved_bounds(y0=ctrl[6:1], y1=trt[6:1], D_ot = 4:6, D_ct = 4:6, tol = 1e-8)

# independence bounds
improved_bounds(y0=ctrl[6:1], y1=trt[6:1], D_ot = 1:6, D_ct = 1:6, tol = 1e-8)

