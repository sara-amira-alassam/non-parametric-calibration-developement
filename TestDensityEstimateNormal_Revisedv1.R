# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(14)

library(carbondate)

###############################################################################
# Create some observed data from the clusters according to probabilities

num_observations <- 100
c14_sig <- rep(25, num_observations)

# Choose base function of calendar age
weights_true <- c(0.1, 0.4, 0.5)
weights_true <- weights_true / sum(weights_true)
phi_true <- c(3500, 4200, 5000)
tau_true <- 1 / c(200, 100, 300)^2
n_clust_true <- length(weights_true)

true_density = data.frame(x = intcal20$calendar_age, y = 0)
true_density$y <- 0
for(i in 1:length(weights_true)) {
  true_density$y <- true_density$y +
    weights_true[i] * dnorm(
      true_density$x, mean = phi_true[i], sd = 1/sqrt(tau_true[i]))
}

cluster_identifiers_true <- sample(
  1:n_clust_true, num_observations, replace = TRUE, prob = weights_true)
calendar_ages_true <- rnorm(
  num_observations,
  mean = phi_true[cluster_identifiers_true],
  sd = 1 / sqrt(tau_true[cluster_identifiers_true]))
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
  n_iter = 10000,
  n_thin = 50)

###############################################################################
# Implement the Walker version of the DPMM
walker_output <- WalkerBivarDirichlet(
  c14_determinations = c14_ages,
  c14_sigmas = c14_sig,
  calibration_curve=intcal20,
  n_iter = 10000,
  n_thin = 50)

###############################################################################
# Plot results
# Create a layout with 2/3 showing the predictive density 1/3*1/2 showing the
# number of clusters for each method
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(3, 3), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(
  output_data = list(walker_output, polya_urn_output),
  true_density = true_density,
  n_posterior_samples = 5000)

PlotNumberOfClusters(polya_urn_output)

PlotNumberOfClusters(walker_output)

par(mfrow = c(1,1))
