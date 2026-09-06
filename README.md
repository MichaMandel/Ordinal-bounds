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
- Demonstrates results on simulated and real data

## Files

- generate colored matrix w title.R  
  Scripts for generating figures using colored matrices

- plot_colored_matrix.R  
  Function for plotting colored K×K matrices

- latex_joint_independence_table.R  
  Generates LaTeX tables for joint and independence distributions

- ordinal_effect_measures.R  
  Computes η, τ and their independence counterparts

- tables for ordinal bounds.R  
  Example matrices and numerical illustrations

- stroke data analysis.R  
  Analysis of stroke trial data (modified Rankin scale)

## Usage

### Compute causal estimands

```r
source("ordinal_effect_measures.R")

mat <- matrix(
  c(0.02, 0.04, 0.06,
    0.09, 0.07, 0.48,
    0.03, 0.02, 0.19),
  nrow = 3, byrow = TRUE
)

ordinal_effect_measures(mat)
