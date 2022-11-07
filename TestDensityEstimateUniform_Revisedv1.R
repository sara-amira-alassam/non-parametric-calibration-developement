# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(83)

library(carbondate)
###############################################################################
# Create some observed data from a uniform distribution

num_observations <- 100
c14_sig <- rep(25, num_observations)

# Create uniform distribution
begin = 6000
end = 6500
true_density = data.frame(
  x = intcal20$calendar_age,
  y = dunif(x = intcal20$calendar_age, min = begin, max = end))
calendar_ages_true <- runif(num_observations, min = begin, max = end)
hist(calendar_ages_true, breaks = 20)

# Create some radiocarbon determinations
interpolated_calibration_curve = InterpolateCalibrationCurve(
  new_calendar_ages = calendar_ages_true, calibration_curve = intcal20)
interpolated_c14_age <- interpolated_calibration_curve$c14_age
interpolated_c14_sig <- interpolated_calibration_curve$c14_sig

# Sample some calibration curve values
xcalcurve <- rnorm(num_observations, interpolated_c14_age, interpolated_c14_sig)
c14_ages <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)

###############################################################################
# Implement the Neal version of the DPMM
polya_urn_output <- PolyaUrnBivarDirichlet(
  c14_determinations = c14_ages,
  c14_sigmas = c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 5)

###############################################################################
# Implement the Walker version of the DPMM
walker_output <- WalkerBivarDirichlet(
  c14_determinations = c14_ages,
  c14_sigmas = c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e5,
  n_thin = 5)

###############################################################################
# Plot results
# Create a layout with 2/3 showing the predictive density 1/3*1/2 showing the
# number of clusters for each method
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(3, 3), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(
  output_data = list(walker_output, polya_urn_output),
  n_posterior_samples = 5000,
  show_confidence_intervals = FALSE,
  true_density=true_density)

PlotNumberOfClusters(polya_urn_output)

PlotNumberOfClusters(walker_output)

par(mfrow = c(1,1))
