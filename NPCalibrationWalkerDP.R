# 13th April 2022
# Code to implement NP Bayesian Calibration and Summaristion of Related 14C Determinations
set.seed(8)

library(carbondate)

###############################################################################
# Set parameters - Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(
  CalibrateSingleDetermination,
  kerr$c14_ages,
  kerr$c14_sig,
  MoreArgs = list(calibration_curve=intcal20))
inittheta <- intcal20$calendar_age[apply(initprobs, 2, which.max)]

maxrange <- max(inittheta) - min(inittheta)

# Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
# E[tau] = (1/100)^2 Var[tau] = (1/100)^4
# Interval for sigma2 is approx 1/c(nu2/nu1-2*nu2^2/nu1, nu2/nu1+2*nu2^2/nu1)
tempspread <- 0.1 * mad(inittheta)
tempprec <- 1 / (tempspread)^2
nu1 <- 0.25
nu2 <- nu1 / tempprec

lambda <- (100 / maxrange)^2 # Each muclust ~ N(mutheta, sigma2/lambda)

###############################################################################
# Perform the MCMC update

walker_temp <- WalkerBivarDirichlet(
  c14_determinations = kerr$c14_ages,
  c14_uncertainties = kerr$c14_sig,
  calibration_curve=intcal20,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2,
  alpha_shape = 1,
  alpha_rate = 1,
  n_iter = 1e5,
  n_thin = 10,
  slice_width = max(1000, diff(range(kerr$c14_ages)) / 2),
  slice_multiplier = 10,
  n_clust = 10)

###############################################################################
# Plot results

# Create a layout with 2/3 showing the predictive density 1/3 showing the number
# of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix, heights = c(1), widths = c(10, 4.5))

PlotCalendarAgeDensity(
  c14_determinations = kerr$c14_ages,
  c14_uncertainties = kerr$c14_sig,
  calibration_curve = intcal20,
  output_data = walker_temp,
  n_posterior_samples = 5000)

PlotNumberOfClusters(output_data = walker_temp)

# New plot for a single determination
par(mfrow = c(1,1))

PlotIndividualCalendarAgeDensity(
  ident=10,
  c14_determinations = kerr$c14_ages,
  c14_uncertainties = kerr$c14_sig,
  calibration_curve = intcal20,
  output_data = walker_temp)
