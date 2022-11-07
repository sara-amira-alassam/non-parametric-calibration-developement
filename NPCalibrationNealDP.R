# 13th April 2022 - Edited October 2022
# Code to implement NP Bayesian Calibration and Summarisation of Related 14C
# Determinations

# NOTE:
# This version uses the Polya Urn method of sampling form the DP
# It is slower so I suggest using the Walker version

set.seed(8)

library(carbondate)

###############################################################################
# Perform the MCMC update
polya_urn_output <- PolyaUrnBivarDirichlet(
  c14_determinations = kerr$c14_ages,
  c14_sigmas = kerr$c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 50)

###############################################################################
# Plot results

# Create a layout with 2/3 showing the predictive density 1/3 showing the number
# of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix, heights = c(1), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(polya_urn_output, n_posterior_samples = 5000)

PlotNumberOfClusters(polya_urn_output)

# New plot for a single determination
par(mfrow = c(1,1))

PlotCalendarAgeDensityIndividualSample(15, polya_urn_output)
