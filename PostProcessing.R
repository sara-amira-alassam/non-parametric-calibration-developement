# This code will plot the predictive density for Walker

plot_final_graphs <- function(Temp, npostsum, calcurve, lambda, nu1, nu2, postden, postdenCI, SPD, x, xsig) {
  # Find predictive density
  npost <- dim(Temp$delta)[1]
  nburn <- floor(npost / 2)
  SPDcol <- grey(0.1, alpha = 0.5)

  # Create the range over which to plot the density
  # Note since tempx is not seq(T_a, T_b, 1) (i.e. spaced every year) then
  # sum(postdenCI) != 1
  # We should have sum(postdenCI) * gridgap = 1 (where gridgap is the spacing in tempx)
  tempx <- seq(
    floor(min(Temp$theta, na.rm = TRUE)),
    ceiling(max(Temp$theta, na.rm = TRUE)),
    by = 1,
  )

  # Now choose the ids of the posterior sample
  sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost - nburn))

  # Create a matrix where each column is the density for a particular sample id
  # We can then find the mean along each row
  postDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) {
    WalkerFindpred(x,
      w = out$w[[i]], phi = out$phi[[i]], tau = out$tau[[i]],
      muphi = out$muphi[i], lambda = lambda, nu1 = nu1, nu2 = nu2
    )
  },
  out = Temp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2
  )

  # Find CI and mean
  postdenCI <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
  postden <- apply(postDmat, 1, mean)

  # Create a layout which has 2/3 of the plot showing the predictive density and 1/3 showing the number of clusters
  layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  layout(
    mat = layout.matrix,
    heights = c(1), # Heights of the two rows
    widths = c(10, 4.5)
  ) # Widths of the two columns

  # Plot the predictive joint density
  xlim <- rev(range(tempx))
  ylim <- range(x) + c(-2, 2) * quantile(xsig, 0.9)
  denscale <- 1.2
  par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
  plot(
    calcurve$calage,
    calcurve$c14age,
    col = "blue",
    ylim = ylim,
    xlim = xlim,
    xlab = "Calendar Age (cal yr BP)",
    ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
    type = "l",
    main = expression(paste(""^14, "C Calibration - Walker DP")),
  )
  calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
  calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
  lines(calcurve$calage, calcurve$ub, lty = 2, col = "blue")
  lines(calcurve$calage, calcurve$lb, lty = 2, col = "blue")
  polygon(
    c(rev(calcurve$calage), calcurve$calage),
    c(rev(calcurve$lb), calcurve$ub),
    col = rgb(0, 0, 1, .3),
    border = NA,
  )
  rug(x, side = 2)

  # Plot the SPD and DPMM density along the bottom
  par(new = TRUE, las = 1)
  plot(tempx, postden,
       lty = 1, col = "purple", type = "n",
       ylim = c(0, denscale * max(postdenCI)), xlim = xlim,
       axes = FALSE, xlab = NA, ylab = NA
  )
  polygon(c(SPD$calage, rev(SPD$calage)), c(SPD$prob, rep(0, length(SPD$prob))),
          border = NA, col = SPDcol
  )
  lines(tempx, postden, col = "purple")
  lines(tempx, postdenCI[1, ], col = "purple", lty = 2)
  lines(tempx, postdenCI[2, ], col = "purple", lty = 2)
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, line = -1.1)

  legend("topright",
    legend = c("IntCal20", "NP Density Estimate - Walker", "95% Prob. Interval", "SPD Estimate"),
    lty = c(1, 1, 2, -1), pch = c(NA, NA, NA, 15),
    col = c("blue", "purple", "black", SPDcol), cex = 0.8, pt.cex = 2
  )

  # Plot the number of clusters
  WalkerNClust <- apply(Temp$delta, 1, function(x) length(unique(x)))
  WalkerNClust <- WalkerNClust[nburn:npost]
  hist(WalkerNClust,
    xlab = "Number of Clusters", main = "",
    probability = TRUE, breaks = seq(0.5, max(WalkerNClust) + 1, by = 1)
  )
  mtext(paste0("(", letters[2], ")"),
    side = 3, adj = 1,
    line = -1.1
  )
}
