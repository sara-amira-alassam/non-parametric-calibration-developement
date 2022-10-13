# 13th April 2022
# Code to implement NP Bayesian Calibration and Summaristion of Related 14C Determinations
set.seed(14)

# Read in the necessary functions
source("WalkerDirichletMixtureUpdateFunsFinal.R") # This also reads in the slice sampling SliceUpdateFuns.R
source("NealDirichletMixtureMasterFunctionsFinal.R")
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")
source("PostProcessing.R")

library(carbondate)

# Read in data
# x - c14ages
# xsig - corresponding 1 sigma
Kerr <- read.csv("RealDatasets/kerr2014sss_sup.csv", header = FALSE, sep = ",")
x <- Kerr[, 3]
xsig <- Kerr[, 4]

#############################################################################
# Now choose hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(cprshape, cprrate)
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- 1
cprrate <- 1

#### Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(
  carbondate::CalibrateSingleDetermination,
  x,
  xsig,
  MoreArgs = list(calibration_curve=carbondate::intcal20))
inittheta <- intcal20$calendar_age[apply(initprobs, 2, which.max)]
# Choose A and B from range of theta
A <- median(inittheta)
B <- 1 / (max(inittheta) - min(inittheta))^2
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


# Choose number of iterations for sampler
niter <- 1000
nthin <- 5 # Don't choose too high, after burn-in we have (niter/nthin)/2 samples from posterior to potentially use
npostsum <- 5000 # Current number of samples it will draw from this posterior to estimate fhat (possibly repeats)

WalkerTemp <- carbondate::WalkerBivarDirichlet(
  c14_determinations = x,
  c14_uncertainties = xsig,
  calibration_curve=intcal20,
  lambda = lambda,
  nu1 = nu1,
  nu2 = nu2,
  A = A,
  B = B,
  cprshape = cprshape,
  cprrate = cprrate,
  n_iter = niter,
  n_thin = nthin,
  calendar_ages = inittheta,
  slice_width = max(1000, diff(range(x)) / 2),
  slice_multiplier = 10,
  kstar = 10)

SPD <- carbondate::FindSPD(
  calendar_age_range=floor(range(WalkerTemp$calendar_ages)),
  c14_determinations=x,
  c14_uncertainties=xsig,
  calibration_curve=carbondate::intcal20)

post_process_and_plot(WalkerTemp, NULL, SPD, NULL, npostsum, intcal20, lambda, nu1, nu2, x, xsig)

# If we want to plot e.g. the posterior calendar age density against the curve then we can run the below
# ident is the determination you want to calibrate
plotindpost(WalkerTemp, ident = 10, y = x, er = xsig, calibration_curve = intcal20)
