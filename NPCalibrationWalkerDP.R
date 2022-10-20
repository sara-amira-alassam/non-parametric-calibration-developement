# 13th April 2022
# Code to implement NP Bayesian Calibration and Summaristion of Related 14C Determinations
set.seed(8)

library(carbondate)

# Read in data
Kerr <- read.csv("RealDatasets/kerr2014sss_sup.csv", header = FALSE, sep = ",")
c14_ages <- Kerr[, 3]
c14_sig <- Kerr[, 4] # corresponding 1 sigma

#### Updated adaptive version
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
# Interval for sigma2 is approx 1/ c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 + 2*nu2^2/nu1)
tempspread <- 0.1 * mad(inittheta)
tempprec <- 1 / (tempspread)^2
nu1 <- 0.25
nu2 <- nu1 / tempprec

# Setup the NP method
lambda <- (100 / maxrange)^2 # Each muclust ~ N(mutheta, sigma2/lambda)

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

carbondate::PlotCalendarAgeDensity(
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  output_data = walker_temp,
  n_posterior_samples = 5000,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2,
)

carbondate::PlotNumberOfClusters(output_data = walker_temp)

carbondate::PlotIndividualCalendarAgeDensity(
  ident=10,
  c14_determinations = c14_ages,
  c14_uncertainties = c14_sig,
  calibration_curve = intcal20,
  output_data = walker_temp)
