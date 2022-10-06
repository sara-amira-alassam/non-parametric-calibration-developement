
find_spd_estimate <- function(yrange, x, xsig, calcurve) {
  yfromto <- seq(max(0, yrange[1] - 400), min(50000, yrange[2] + 400), by = 1)

  # Find the calibration curve mean and sd over the yrange
  CurveR <- FindCalCurve(yfromto, calcurve)

  # Now we want to apply to each radiocarbon determination
  # Matrix where each column represents the posterior probability of each theta in yfromto
  # indprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = CurveR$curvemean, calsig = CurveR$curvesd))
  calibration_data = data.frame(c14_age = CurveR$curvemean, c14_sig = CurveR$curvesd)
  indprobs <- mapply(carbondate::CalibrateSingleDetermination, x, xsig, MoreArgs = list(calibration_data = calibration_data))

  # Find the SPD estimate (save as dataframe)
  SPD <- data.frame(
    calage = yfromto,
    prob = apply(indprobs, 1, sum) / dim(indprobs)[2]
  )

  return(SPD)
}


post_process_and_plot <- function(
    WalkerTemp, NealTemp, SPD, true_density, npostsum, calcurve, lambda, nu1, nu2, x, xsig,
    xlimscal = 1, ylimscal=1, denscale=3) {
  SPD_colour <- grey(0.1, alpha = 0.5)
  calibration_curve_colour <- "blue"
  walker_colour <- "purple"
  neal_colour <- "forestgreen"
  true_density_colour <- "red"

  tempx <- create_range_to_plot_density(WalkerTemp, NealTemp)

  xlim <- scale_limit(rev(range(tempx)), xlimscal)
  ylim <- scale_limit(range(x) + c(-2, 2) * quantile(xsig, 0.9), ylimscal)

  create_plot_layout(WalkerTemp, NealTemp)
  plot_calibration_curve_and_data(xlim, ylim, calcurve, x, calibration_curve_colour)

  # Plot the SPD and DPMM density along the bottom
  plot_SPD_estimate_on_current_plot(SPD, SPD_colour, denscale, xlim)
  if (!is.null(WalkerTemp)) {
    add_walker_density_estimate(walker_colour, WalkerTemp, tempx, npostsum, lambda, nu1, nu2)
  }
  if (!is.null(NealTemp)) {
    add_neal_density_estimate(neal_colour, NealTemp, tempx, npostsum, lambda, nu1, nu2)
  }
  if (!is.null(true_density)) {
    lines(true_density$x, true_density$y, col = true_density_colour)
  }
  add_legend_to_density_plot(
    WalkerTemp, NealTemp, true_density, calibration_curve_colour, walker_colour, neal_colour,
    true_density_colour, SPD_colour)
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, line = -1.1)

  plot_index <- 2
  if (!is.null(WalkerTemp)) {
    plot_number_of_walker_clusters(WalkerTemp)
    plot_index <-3
  }
  if (!is.null(NealTemp)) {
    plot_number_of_neal_clusters(NealTemp, plot_index)
  }
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


scale_limit <- function(lim, limscal) {
  lim <- lim + c(1, -1) * diff(lim) * (1 - limscal)
  return(lim)
}


create_plot_layout <- function(WalkerTemp, NealTemp) {
  if (!is.null(WalkerTemp) && !is.null(NealTemp)) {
    # Create a layout with 2/3 showing the predictive density 1/3*1/2 showing the number of clusters
    # for each method
    layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)
    layout(
      mat = layout.matrix,
      heights = c(3, 3),
      widths = c(10, 4.5)
    )
  } else {
    # Create a layout with 2/3 showing the predictive density 1/3 showing the number of clusters
    layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
    layout(
      mat = layout.matrix,
      heights = c(1),
      widths = c(10, 4.5)
    )
  }
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
    main = expression(paste(""^14, "C Calibration")),
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


add_neal_density_estimate <- function(neal_colour, NealTemp, tempx, npostsum, lambda, nu1, nu2) {
  postDmat <- find_density_per_sample_id_neal(NealTemp, tempx, npostsum, lambda, nu1, nu2)
  # Find CI and mean of density along each row
  postdenCI_neal <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
  postden_neal <- apply(postDmat, 1, mean)

  lines(tempx, postden_neal, col = neal_colour)
  lines(tempx, postdenCI_neal[1, ], col = neal_colour, lty = 2)
  lines(tempx, postdenCI_neal[2, ], col = neal_colour, lty = 2)
}


find_density_per_sample_id_neal <- function(NealTemp, tempx, npostsum, lambda, nu1, nu2) {
  npost <- dim(NealTemp$theta)[1]
  nburn <- floor(npost / 2)

  # Now choose the ids of the posterior sample
  sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost - nburn))

  # Create a matrix where each column is the density for a particular sample id
  postDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) {
    NealFindpred(x,
                 c = out$c[i, ], phi = out$phi[[i]], tau = out$tau[[i]],
                 alpha = out$alpha[i], muphi = out$muphi[i],
                 lambda = lambda, nu1 = nu1, nu2 = nu2
    )
  },
  out = NealTemp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2
  )
  return(postDmat)
}


add_legend_to_density_plot <- function(
    WalkerTemp, NealTemp, true_density, calibration_curve_colour, walker_colour, neal_colour,
    true_density_colour, SPD_colour) {
  legend_labels = "IntCal20"
  lty = 1
  pch = NA
  col = calibration_curve_colour
  if (!is.null(true_density)) {
    legend_labels <- c(legend_labels, "True density")
    lty <- c(lty, 1)
    pch <- c(pch, NA)
    col <- c(col, true_density_colour)
  }
  if (!is.null(WalkerTemp)) {
    legend_labels <- c(legend_labels, "Walker DP", "Walker 95% prob interval")
    lty <- c(lty, 1, 2)
    pch <- c(pch, NA, NA)
    col <- c(col, walker_colour, walker_colour)
  }
  if (!is.null(NealTemp)) {
    legend_labels <- c(legend_labels, "Neal DP", "Neal 95% prob interval")
    lty <- c(lty, 1, 2)
    pch <- c(pch, NA, NA)
    col <- c(col, neal_colour, neal_colour)
  }

  legend_labels <- c(legend_labels, "SPD Estimate")
  lty <- c(lty, -1)
  pch <- c(pch, 15)
  col <- c(col, SPD_colour)


  legend("topright", legend = legend_labels, lty = lty, pch = pch, col = col)
}


plot_number_of_walker_clusters <- function(WalkerTemp) {
  npost <- dim(WalkerTemp$delta)[1]
  nburn <- floor(npost / 2)

  WalkerNClust <- apply(WalkerTemp$delta, 1, function(x) length(unique(x)))
  WalkerNClust <- WalkerNClust[nburn:npost]
  hist(WalkerNClust,
       xlab = "Number of Clusters", main = "Walker - Slice Sample DP",
       probability = TRUE, breaks = seq(0.5, max(WalkerNClust) + 1, by = 1)
  )
  mtext(paste0("(", letters[2], ")"),
        side = 3, adj = 1,
        line = -1.1
  )
}


plot_number_of_neal_clusters <- function(NealTemp, plot_index) {
  npost <- dim(NealTemp$theta)[1]
  nburn <- floor(npost / 2)

  NealNClust <- apply(NealTemp$c, 1, max)
  NealNClust <- NealNClust[nburn:npost]
  hist(NealNClust,
       xlab = "Number of Clusters", main = "Neal - P\u{F2}lya Urn DP",
       probability = TRUE, breaks = seq(0.5, max(NealNClust) + 1, by = 1)
  )
  mtext(paste0("(", letters[plot_index], ")"),
        side = 3, adj = 1,
        line = -1.1
  )
}
