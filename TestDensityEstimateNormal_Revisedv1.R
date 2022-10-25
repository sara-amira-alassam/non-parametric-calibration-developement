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
# Set parameters - Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(
  CalibrateSingleDetermination,
  c14_ages,
  c14_sig,
  MoreArgs = list(calibration_curve = intcal20))
inittheta <- intcal20$calendar_age[apply(initprobs, 2, which.max)]

maxrange <- max(inittheta) - min(inittheta)

# Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
# E[tau] = (1/100)^2 Var[tau] = (1/100)^4
# Interval for sigma2 is approx 1/c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 +2*nu2^2/nu1)
tempspread <- 0.1 * mad(inittheta)
tempprec <- 1 / (tempspread)^2
nu1 <- 0.25
nu2 <- nu1 / tempprec
lambda <- (100 / maxrange)^2 # Each muclust ~ N(mutheta, sigma2/lambda)


###############################################################################
# Implement the Neal version of the DPMM
neal_temp <- BivarGibbsDirichletwithSlice(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2,
  alpha_shape = 1,
  alpha_rate = 1,
  n_iter = 1000,
  n_thin = 5,
  slice_width = max(1000, diff(range(c14_ages)) / 2),
  slice_multiplier = 10,
  n_clust = 10)

###############################################################################
# Implement the Walker version of the DPMM
walker_temp <- WalkerBivarDirichlet(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve=intcal20,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2,
  alpha_shape = 1,
  alpha_rate = 1,
  n_iter = 1000,
  n_thin = 5,
  slice_width = max(1000, diff(range(c14_ages)) / 2),
  slice_multiplier = 10,
  n_clust = 10)


###############################################################################
# Plot results
# Create a layout with 2/3 showing the predictive density 1/3*1/2 showing the
# number of clusters for each method
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(3, 3), widths = c(10, 4.5))

PlotCalendarAgeDensity(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  output_data = list(walker_temp, neal_temp),
  true_density = true_density,
  n_posterior_samples = 500,
  show_confidence_intervals = FALSE)

PlotNumberOfClusters(output_data = neal_temp)

PlotNumberOfClusters(output_data = walker_temp)

par(mfrow = c(1,1))
