# AEPPMx
A covariate-dependent product partition model (PPMx) for aggregate adverse event (AE) monitoring across clinical trials

# Introduction
This R package implements a Bayesian nonparametric model for aggregate adverse event (AE) monitoring across clinical trials. The method leverages a covariate-dependent product partition model (PPMx) to integrate information from both historical and ongoing studies, accommodating data at varying levels of granularity—from individual-level covariates to study-level summaries. A novel pairwise similarity measure is used to guide the clustering of experimental units with similar covariate profiles, enhancing the precision of AE rate estimation. The framework supports blinded real-time safety monitoring with seamless transition to unblinded analyses when needed. Case studies and simulation tools are included to demonstrate safety signal detection and risk assessment under diverse trial scenarios.

# R Package tutorial
See the two example files:

1. **example_blind.R**

2. **example_unblind.R**
