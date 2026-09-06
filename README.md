# Ordinal Causal Bounds

This repository contains R code accompanying a research project on causal inference with ordinal outcomes, focusing on bounding causal estimands and studying dependence conditions under which tighter bounds are valid.

## Overview

For ordinal outcomes, standard causal estimands such as the average treatment effect are not well suited. Instead, we consider:

- η = P(Y(1) > Y(0))
- τ = P(Y(1) ≥ Y(0))

These quantities are generally not identifiable, but can be bounded using the marginal distributions of the potential outcomes.

This project:

- Computes nonparametric bounds
- Computes independence-based bounds
- Introduces and studies diagonal tail dominance (DTD)
- Provides visualization tools using colored matrices
- Demonstrates results on theoretical and real data

## Files

- analyze_ordinal_joint.R
  
  Calculates functional from a joint probability matrix
  
- analyze_ordinal_marginals.R
  
  Computes bounds from marginal distributions
  
- generate_colored_matrix_w_title.R
  
  Generates colored matrix of different estimands and bounds
  
- joint_gaussian_copula.R
  
  Generates an ordinal gussian copula for given marginals and \rho
  
- latex_joint_independence_table.R
  
  Generates LaTeX tables for joint and independence distributions

- stroke data analysis.R
  
  Analysis of stroke trial data (modified Rankin scale)

- tables for ordinal bounds.R
  
  Example matrices and numerical illustrations


