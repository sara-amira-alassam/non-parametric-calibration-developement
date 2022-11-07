# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(17)

library(carbondate)

###############################################################################
# Read in population data and simulate some radiocarbon determinations from it

# Read in Mexico population data
Mexico <- read.csv("ContrerasMeadows/MexicoPop.csv", header = TRUE)
plot(Mexico$calBP, Mexico$Pop, type = "l")
mintheta <- ceiling(min(Mexico$calBP))
maxtheta <- floor(max(Mexico$calBP))
popinterp <- approx(x = Mexico$calBP, y = Mexico$Pop, xout = mintheta:maxtheta)
# Re-scale for density and convert to data frame for later plotting
popinterp <- data.frame(x = popinterp$x, y = popinterp$y / sum(popinterp$y))

# Create simulated calendar ages
num_observations <- 500
true_calendar_ages <- sample(
  popinterp$x, num_observations, replace = TRUE, prob = popinterp$y)
hist(true_calendar_ages, breaks = 30)

# Create some radiocarbon determinations x
interpolated_calibration_curve = InterpolateCalibrationCurve(
  new_calendar_ages = true_calendar_ages, calibration_curve = intcal20)
interpolated_c14_age <- interpolated_calibration_curve$c14_age
interpolated_c14_sig <- interpolated_calibration_curve$c14_sig

# Sample some calibration curve values
xcalcurve <- rnorm(num_observations, interpolated_c14_age, interpolated_c14_sig)
c14_sig <- round(runif(num_observations, min = 20, max = 40), 0)
c14_ages <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)

###############################################################################
# Implement the Polya Urn version of the DPMM
polya_urn_output <- PolyaUrnBivarDirichlet(
  c14_determinations = c14_ages,
  c14_sigmas = c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 50)

###############################################################################
# Implement the Walker version of the DPMM
walker_output <- WalkerBivarDirichlet(
  c14_determinations = c14_ages,
  c14_sigmas = c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e5,
  n_thin = 50)

###############################################################################
# Plot results
# Create a layout with 2/3 showing the predictive density 1/3*1/2 showing the
# number of clusters for each method
layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(3, 3), widths = c(10, 4.5))

PlotPredictiveCalendarAgeDensity(
  list(walker_output, polya_urn_output),
  n_posterior_samples = 5000,
  true_density = popinterp)

PlotNumberOfClusters(polya_urn_output)

PlotNumberOfClusters(walker_output)

par(mfrow = c(1,1))
