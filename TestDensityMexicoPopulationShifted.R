# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(17) # Same seed as non-shifted so should generate same cal age sample (just shifted)

# Read in the necessary functions
source("WalkerDirichletMixtureUpdateFunsFinal.R") # This also reads in the slice sampling SliceUpdateFuns.R
source("NealDirichletMixtureMasterFunctionsFinal.R")
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")
source("PostProcessing.R")

# Read in IntCal20 curve
calcurve <- read.table("Curves/intcal20.14c", sep = ",", header = FALSE, skip = 11)
names(calcurve) <- c("calage", "c14age", "c14sig", "Delta14C", "DeltaSigma")


# Read in Mexico populatin data
Mexico <- read.csv("ContrerasMeadows/MexicoPop.csv", header = TRUE)
Mexico$calBP <- Mexico$calBP + 6500
plot(Mexico$calBP, Mexico$Pop, type = "l")
mintheta <- ceiling(min(Mexico$calBP))
maxtheta <- floor(max(Mexico$calBP))
popinterp <- approx(x = Mexico$calBP, y = Mexico$Pop, xout = mintheta:maxtheta)
popinterp$y <- popinterp$y / sum(popinterp$y) # Rescale for density
# Create simulated thetas
nobs <- 500 # Number of observations
thetatrue <- sample(popinterp$x, nobs, replace = TRUE, prob = popinterp$y)
xsig <- round(runif(nobs, min = 20, max = 40), 0)
hist(thetatrue, breaks = 30)

#### Now create some radiocarbon determinations x

# Interpolate calobration curve mean and sd at theta values
calinterp <- FindCal(thetatrue, calmu = calcurve$c14age, caltheta = calcurve$calage, calsig = calcurve$c14sig)

# Sample some calibration curve values
xcalcurve <- rnorm(nobs, calinterp$mu, calinterp$sigma)

x <- rnorm(nobs, mean = xcalcurve, sd = xsig)

###############################################
# Now run the Walker Sampler
###############################################
# Choose number of iterations for sampler
niter <- 1000
nthin <- 5
npostsum <- 5000

############################################################################
# Now choose fixed DP hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(alphaprshape, alphaprrate)
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- alphaprshape <- 1
cprrate <- alphaprrate <- 1

#### Updated adaptive version
# Prior on mu theta for DP - very uninformative based on observed data
initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = calcurve$c14age, calsig = calcurve$c14sig))
inittheta <- calcurve$calage[apply(initprobs, 2, which.max)]
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

#####################################################################
# Implement the Neal version of the DPMM
NealTemp <- BivarGibbsDirichletwithSlice(
  x = x, xsig = xsig,
  lambda = lambda, nu1 = nu1, nu2 = nu2,
  A = A, B = B,
  mualpha = NA, sigalpha = NA,
  alphaprshape = alphaprshape, alphaprrate = alphaprrate,
  niter = niter, nthin = nthin, theta = inittheta,
  w = max(1000, diff(range(x)) / 2), m = 10,
  calcurve = calcurve, nclusinit = 10
)

#####################################################################
# Implement the Walker version of the DPMM
WalkerTemp <- WalkerBivarDirichlet(
  x = x, xsig = xsig,
  lambda = lambda, nu1 = nu1, nu2 = nu2,
  A = A, B = B,
  cprshape = cprshape, cprrate = cprrate,
  niter = niter, nthin = nthin, theta = inittheta,
  slicew = max(1000, diff(range(x)) / 2), m = 10,
  calcurve = calcurve, kstar = 10
)
Temp <- WalkerTemp


SPD <- find_spd_estimate(yrange=floor(range(WalkerTemp$theta)), x, xsig, calcurve)

post_process_and_plot(
  WalkerTemp, NealTemp, SPD, popinterp, npostsum, calcurve, lambda, nu1, nu2, x, xsig,
  xlimscal=1, ylimscal=1, denscale=6)
