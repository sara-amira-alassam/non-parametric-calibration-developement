# This code will plot the predictive density for Walker

post_process_and_plot <- function(WalkerTemp, NealTemp, SPD, npostsum, calcurve, lambda, nu1, nu2, postden, postdenCI, x, xsig) {
  SPD_colour <- grey(0.1, alpha = 0.5)
  calibration_curve_colour <- "blue"
  walker_colour <- "purple"
  neal_colour <- "orange"

  tempx <- create_range_to_plot_density(WalkerTemp, NealTemp)

  xlim <- rev(range(tempx))
  ylim <- range(x) + c(-2, 2) * quantile(xsig, 0.9)
  denscale <- 2.5

  create_plot_layout()
  plot_calibration_curve_and_data(xlim, ylim, calcurve, x, calibration_curve_colour)

  # Plot the SPD and DPMM density along the bottom
  plot_SPD_estimate_on_current_plot(SPD, SPD_colour, denscale, xlim)
  if (!is.null(WalkerTemp)) {
    add_walker_density_estimate(walker_colour, WalkerTemp, tempx, npostsum, lambda, nu1, nu2)
  }
  add_legend_to_density_plot(WalkerTemp, NealTemp, calibration_curve_colour, walker_colour, neal_colour, SPD_colour)
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, line = -1.1)

  plot_number_of_walker_clusters(WalkerTemp)
}


create_range_to_plot_density <- function(WalkerTemp, NealTemp) {
  # Create the range over which to plot the density
  # Note since tempx is not seq(T_a, T_b, 1) (i.e. spaced every year) then
  # sum(postdenCI) != 1
  # We should have sum(postdenCI) * gridgap = 1 (where gridgap is the spacing in tempx)
  if (!is.null(WalkerTemp)) {
    tempx <- seq(
      floor(min(WalkerTemp$theta, na.rm = TRUE)),
      ceiling(max(WalkerTemp$theta, na.rm = TRUE)),
      by = 1,
    )
  } else {
    tempx <- seq(
      floor(min(NealTemp$theta, na.rm = TRUE)),
      ceiling(max(NealTemp$theta, na.rm = TRUE)),
      by = 1,
    )
  }
  return(tempx)
}


create_plot_layout <- function() {
  # Create a layout with 2/3 showing the predictive density 1/3 showing the number of clusters
  layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
  layout(
    mat = layout.matrix,
    heights = c(1),
    widths = c(10, 4.5)
  )
}


plot_calibration_curve_and_data <- function(xlim, ylim, calcurve, x, calibration_curve_colour) {
  par(mar = c(5, 4.5, 4, 2) + 0.1, las = 1)
  plot(
    calcurve$calage,
    calcurve$c14age,
    col = calibration_curve_colour,
    ylim = ylim,
    xlim = xlim,
    xlab = "Calendar Age (cal yr BP)",
    ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
    type = "l",
    main = expression(paste(""^14, "C Calibration - Walker DP")),
  )
  calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
  calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
  lines(calcurve$calage, calcurve$ub, lty = 2, col = calibration_curve_colour)
  lines(calcurve$calage, calcurve$lb, lty = 2, col = calibration_curve_colour)
  polygon(
    c(rev(calcurve$calage), calcurve$calage),
    c(rev(calcurve$lb), calcurve$ub),
    col = rgb(0, 0, 1, .3),
    border = NA,
  )
  rug(x, side = 2)
}


plot_SPD_estimate_on_current_plot <- function(SPD, SPD_colour, denscale, xlim) {
  par(new = TRUE)
  plot(SPD$calage, SPD$prob,
       lty = 1, col = SPD_colour, type = "n",
       ylim = c(0, denscale * max(SPD$prob)), xlim = xlim,
       axes = FALSE, xlab = NA, ylab = NA
  )
  polygon(c(SPD$calage, rev(SPD$calage)), c(SPD$prob, rep(0, length(SPD$prob))),
          border = NA, col = SPD_colour
  )
}


add_walker_density_estimate <- function(walker_colour, WalkerTemp, tempx, npostsum, lambda, nu1, nu2) {
  postDmat <- find_density_per_sample_id_walker(WalkerTemp, tempx, npostsum, lambda, nu1, nu2)
  # Find CI and mean of density along each row
  postdenCI_walker <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
  postden_walker <- apply(postDmat, 1, mean)

  lines(tempx, postden_walker, col = walker_colour)
  lines(tempx, postdenCI_walker[1, ], col = walker_colour, lty = 2)
  lines(tempx, postdenCI_walker[2, ], col = walker_colour, lty = 2)
}


find_density_per_sample_id_walker <- function(WalkerTemp, tempx, npostsum, lambda, nu1, nu2) {
  npost <- dim(WalkerTemp$delta)[1]
  nburn <- floor(npost / 2)

  # Now choose the ids of the posterior sample
  sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost - nburn))

  # Create a matrix where each column is the density for a particular sample id
  postDmat <- apply(
    as.row(sampid),
    2,
    function(i, out, x, lambda, nu1, nu2) {
      WalkerFindpred(x,
                     w = out$w[[i]], phi = out$phi[[i]], tau = out$tau[[i]],
                     muphi = out$muphi[i], lambda = lambda, nu1 = nu1, nu2 = nu2
      )
    },
    out = WalkerTemp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2
  )
  return(postDmat)
}


add_legend_to_density_plot <- function(WalkerTemp, NealTemp, calibration_curve_colour, walker_colour, neal_colour, SPD_colour) {
  legend_labels = "IntCal20"
  lty = 1
  pch = NA
  col = calibration_curve_colour

  if (!is.null(WalkerTemp)) {
    legend_labels <- c(legend_labels, "Walker DP", "Walker 95% prob interval")
    lty <- c(lty, 1, 2)
    pch <- c(pch, NA, NA)
    col <- c(col, walker_colour, walker_colour)
  }

  legend_labels <- c(legend_labels, "SPD Estimate")
  lty <- c(lty, -1)
  pch <- c(pch, 15)
  col <- c(col, SPD_colour)


  legend("topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}


plot_number_of_walker_clusters <- function(Temp) {

  npost <- dim(Temp$delta)[1]
  nburn <- floor(npost / 2)

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
