# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(17)

library(carbondate)

###############################################################################
# Read in population data and simulate some radiocarbon determinations from it

# Read in Mexico population data
Mexico <- read.csv("ContrerasMeadows/MexicoPop.csv", header = TRUE)
Mexico$calBP <- Mexico$calBP + 6500
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
interpolated_calibration_curve = carbondate::InterpolateCalibrationCurve(
  new_calendar_ages = true_calendar_ages, calibration_curve = intcal20)
interpolated_c14_age <- interpolated_calibration_curve$c14_age
interpolated_c14_sig <- interpolated_calibration_curve$c14_sig

# Sample some calibration curve values
xcalcurve <- rnorm(num_observations, interpolated_c14_age, interpolated_c14_sig)
c14_sig <- round(runif(num_observations, min = 20, max = 40), 0)
c14_ages <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)


###############################################################################
# Set parameters - Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(
  carbondate::CalibrateSingleDetermination,
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
neal_temp <- carbondate::BivarGibbsDirichletwithSlice(
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
walker_temp <- carbondate::WalkerBivarDirichlet(
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

carbondate::PlotCalendarAgeDensity(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  output_data = list(walker_temp, neal_temp),
  n_posterior_samples = 500,
  true_density = popinterp)

carbondate::PlotNumberOfClusters(output_data = neal_temp)

carbondate::PlotNumberOfClusters(output_data = walker_temp)

par(mfrow = c(1,1))
