# 30th March 2022
# This file runs the three real life examples of:
# Kerr - Irish Rath n = 255 obs
# Buchanan - Palaeo-Indian demography n = 628 obs
# Armit - Population Decline in Iron Age n=2021 obs
# using the new IntCal20 curve

library("carbondate")

################################################################################
# Read in the data - use one of `armit`, `buchanan` or `kerr`
example_set <- buchanan

###############################################################################
# Perform the MCMC update

walker_output <- WalkerBivarDirichlet(
  c14_determinations = example_set$c14_ages,
  c14_sigmas = example_set$c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 50)

###############################################################################
# Plot results

# Create a layout with 2/3 showing the predictive density 1/3 showing the number
# of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix, heights = c(1), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(walker_output, n_posterior_samples = 5000)

PlotNumberOfClusters(walker_output)

par(mfrow = c(1,1))
