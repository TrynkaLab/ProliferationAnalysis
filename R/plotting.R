#-------------------------------------------------------------------------------
#' Plot initial parameter estimates
#'
#' @param cur.hist histogram object returned by \code{hist()}.
#' @param y.smth numeric vector of smoothed counts (same length as histogram midpoints).
#' @param peak.stats data frame of initial peak estimates from \code{find_initial_peaks}.
#' @param peak.0.lower.bound numeric; lower bound used to locate generation 0.
#' @param peak.x.model logical; whether a peak X is being modelled.
#' @param peak.x.mode numeric; estimated log10 position of peak X, or NULL.
#' @param peak.x.upper.bound numeric; hard upper limit for peak X position, or NULL.
#' @returns NULL invisibly; called for its side effect of producing a plot.
opt_plot_estimates <- function(cur.hist, y.smth, peak.stats, peak.0.lower.bound, peak.x.model, peak.x.mode, peak.x.upper.bound) {

  plot(cur.hist$mids,
       cur.hist$counts,
       main="Initial peak estimates",
       ylab="Count",
       xlab="Bin average log10(CTV)",
       bty="n",
       pch=20,
       col="grey")

  lines(cur.hist$mids, y.smth, new=F, lwd=3, col="blue")
  abline(v=peak.stats[1, 1], col="red", lwd=2)
  abline(v=peak.0.lower.bound, col="red", lty=2)

  for (i in 2:nrow(peak.stats)) {
    abline(v=peak.stats[i, 1], lwd=2, col="grey")
  }

  if (peak.x.model & !is.null(peak.x.mode)) {
    abline(v=peak.x.mode, col="orange", lwd=2)
    abline(v=peak.x.upper.bound, col="orange", lty=2)

    legend("topleft",
           legend=c("Smoothed data", "Gen0", "GenX"),
           fill=c("blue", "red", "orange"),
           bty="n")
  } else {
    legend("topleft",
           legend=c("Smoothed data", "Gen0"),
           fill=c("blue", "red"),
           bty="n")
  }
}

#-------------------------------------------------------------------------------
#' Plot optimization of parameters
#'
#' @param x.mids numeric vector of histogram midpoints (log10 scale).
#' @param y numeric vector of observed counts matching x.mids in length.
#' @param y.pred.optim numeric vector of model-predicted counts matching x.mids.
#' @param main character string plot title (default "Optimal fit vs data").
#' @returns NULL invisibly; called for its side effect of producing a plot.
opt_plot_final <- function(x.mids, y, y.pred.optim, main="Optimal fit vs data") {

  plot(x.mids, y,
       main=main,
       xlab="Bin average log10(CTV)",
       ylab="Count",
       pch=20,
       col="grey", bty="n")
  lines(y.pred.optim ~ x.mids, lwd=2, col="red")

  legend("topleft",
         legend=c("Data", "Optimized fit"),
         fill=c("grey", "red"),
         bty="n")
}

#-------------------------------------------------------------------------------
#' Plot optimization of parameters per peak
#'
#' @param x.mids numeric vector of histogram midpoints (log10 scale).
#' @param y numeric vector of observed counts matching x.mids in length.
#' @param ps data frame of peak statistics as returned by \code{fit_peaks},
#'   containing columns \code{opt_mean}, \code{opt_peak_sd}, \code{peak_area_prop},
#'   and \code{peak_events}.
#' @param main character string plot title (default "Optimal fit vs data").
#' @param peak.x.model logical; if TRUE the last row of \code{ps} is peak X and
#'   its individual density is drawn in orange (default FALSE).
#' @param peak.0.lower.bound numeric; if supplied, a grey dashed vertical line is
#'   drawn at this log10 position marking the generation-0 lower bound (default NULL).
#' @param peak.x.upper.bound numeric; if supplied, a grey dashed vertical line is
#'   drawn at this log10 position marking the hard upper limit for peak X (default NULL).
#' @returns NULL invisibly; called for its side effect of producing a plot.
opt_plot_final_pp <- function(x.mids, y, ps, main="Optimal fit vs data",
                              peak.x.model=FALSE,
                              peak.0.lower.bound=NULL,
                              peak.x.upper.bound=NULL) {

  plot(x.mids, y,
       main=main,
       xlab="Bin average log10(CTV)",
       ylab="Count",
       pch=20,
       col="grey", bty="n")

  binwidth     <- (x.mids[2] - x.mids[1])

  for (i in 1:nrow(ps)) {
    range    <- x.mids
    y.single <- dnorm(range, ps[i,]$opt_mean, ps[i,]$opt_peak_sd)
    y.single <- y.single * ps[i,]$peak_area_prop/100
    y.single <- (y.single*binwidth) * sum(ps$peak_events)

    peak_col <- if (peak.x.model && i == nrow(ps)) "orange" else "black"
    lines(range, y.single, col=peak_col)
  }

  dens <- exp(prolif_model_density(ps$opt_mean, ps$opt_peak_sd, ps$peak_area_prop, x.mids, log=T))
  lines(x.mids, (dens*binwidth) * sum(ps$peak_events), col="red", lwd=2)

  if (!is.null(peak.0.lower.bound)) {
    abline(v=peak.0.lower.bound, col="darkgrey", lty=2)
  }
  if (!is.null(peak.x.upper.bound)) {
    abline(v=peak.x.upper.bound, col="darkgrey", lty=2)
  }

  if (peak.x.model) {
    legend("topleft",
           legend=c("Data", "Optimized fit", "Single peaks", "Peak X"),
           fill=c("grey", "red", "black", "orange"),
           bty="n")
  } else {
    legend("topleft",
           legend=c("Data", "Optimized fit", "Single peaks"),
           fill=c("grey", "red", "black"),
           bty="n")
  }
}

#-------------------------------------------------------------------------------
#' Plot optimization trace showing convergence of parameters
#'
#' @param x numeric vector of x values (histogram midpoints), used to set the
#'   y-axis range of the means panel.
#' @param opt.env environment created by \code{opt_new_env} and populated
#'   during optimisation, containing \code{score}, \code{means}, \code{sd},
#'   and \code{summits} matrices.
#' @returns NULL invisibly; called for its side effect of producing plots.
opt_plot_params <- function(x, opt.env) {

  n     <- ncol(opt.env[["means"]])
  cols  <- palette.colors(n = n, recycle = T)

  #--------------------------------------------------------------------------
  # Score
  plot(1:length(opt.env[["score"]]),
       opt.env[["score"]],
       xlab="Iteration",
       ylab="NLL / RSS",
       pch=20,
       col="grey",
       bty="n",
       main="Neg. Log Lik. / Res. Sum of Sqr.",
       type="l")

  #--------------------------------------------------------------------------
  n.iter <- nrow(opt.env[["means"]])

  # Mean
  plot(1:n.iter,
       opt.env[["means"]][,1, drop=F],
       xlab="Iteration",
       ylab="Mean",
       pch=20,
       col=cols[1],
       bty="n",
       main="Optim - mean",
       type="l",
       ylim=c(min(x), max(x)))

  if (ncol(opt.env[["means"]]) > 1) {
    for (i in 2:ncol(opt.env[["means"]])) {
      lines(1:n.iter,
            opt.env[["means"]][,i,drop=F],
            col=cols[i])
    }
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")

  #--------------------------------------------------------------------------
  # SD
  plot(1:n.iter,
       opt.env[["sd"]][,1],
       xlab="Iteration",
       ylab="Peak SD",
       pch=20,
       col=cols[1],
       bty="n",
       main="Optim - SD",
       type="l",
       ylim=c(min(opt.env[["sd"]]), max(opt.env[["sd"]])))

  if (ncol(opt.env[["sd"]]) > 1) {
    for (i in 2:ncol(opt.env[["sd"]])) {
      lines(1:n.iter,
            opt.env[["sd"]][,i],
            col=cols[i])
    }
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")

  #--------------------------------------------------------------------------
  # Summit
  plot(1:n.iter,
       opt.env[["summits"]][,1],
       xlab="Iteration",
       ylab="Summit",
       pch=20,
       col=cols[1],
       bty="n",
       main="Optim - summit",
       type="l",
       ylim=c(min(opt.env[["summits"]]), max(opt.env[["summits"]])))

  if (ncol(opt.env[["summits"]]) > 1) {
    for (i in 2:ncol(opt.env[["summits"]])) {
      lines(1:n.iter,
            opt.env[["summits"]][,i],
            col=cols[i])
    }
  }

  legend("right",
         legend=paste0("Gen", 0:(n-1)),
         fill=cols,
         bty="n")
}
