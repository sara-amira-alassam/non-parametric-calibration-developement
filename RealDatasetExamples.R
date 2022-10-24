# 30th March 2022
# This file runs the three real life examples of:
# Kerr - Irish Rath n = 255 obs
# Buchanan - Palaeo-Indian demography n = 628 obs
# Armit - Population Decline in Iron Age n=2021 obs
# using the new IntCal20 curve

library("carbondate")

################################################################################
# Read in the data
ExampleSet <- "Armit"

if (ExampleSet == "Kerr") {
  Kerr <- read.csv("RealDatasets/kerr2014sss_sup.csv", header = FALSE, sep = ",")
  c14_ages <- Kerr[, 3]
  c14_sig <- Kerr[, 4] # corresponding 1 sigma
} else if (ExampleSet == "Buchanan") {
  Data <- read.csv("RealDatasets/buchanan2008pde.csv", header = FALSE, sep = ",")
  c14_ages <- Data[, 3]
  c14_sig <- Data[, 4]
} else if (ExampleSet == "Armit") {
  Armit <- read.csv("RealDatasets/armit2014rcc_sd01.csv", header = TRUE, sep = ",")
  # Remove the two observations with missing xsig values
  remove <- which(is.na(Armit$error))
  Armit <- Armit[-remove, ]
  c14_ages <- Armit$X14C.age
  c14_sig <- Armit$error
} else {
  stop("Unknown Example")
}

###############################################################################
# Set parameters - Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(
  carbondate::CalibrateSingleDetermination,
  c14_ages,
  c14_sig,
  MoreArgs = list(calibration_curve=carbondate::intcal20))
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

# Create a layout with 2/3 showing the predictive density 1/3 showing the number
# of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix, heights = c(1), widths = c(10, 4.5))

carbondate::PlotCalendarAgeDensity(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  output_data = walker_temp,
  n_posterior_samples = 5000,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2)

carbondate::PlotNumberOfClusters(output_data = walker_temp)
