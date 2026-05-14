#' @import minpack.lm
NULL

#-------------------------------------------------------------------------------
#' Create a new env for storing optimization trace
#'
#' @param n.peaks integer number of peaks; sets the number of columns in the
#'   logged parameter matrices.
#' @returns an environment with pre-allocated fields: \code{score} (numeric
#'   vector), \code{means}, \code{sd}, and \code{summits} (matrices with
#'   \code{n.peaks} columns).
opt_new_env <- function(n.peaks) {
  opt.env              <- new.env()
  opt.env[["score"]]   <- c()
  opt.env[["means"]]   <- matrix(ncol=n.peaks, nrow=0)
  opt.env[["sd"]]      <- matrix(ncol=n.peaks, nrow=0)
  opt.env[["summits"]] <- matrix(ncol=n.peaks, nrow=0)

  return(opt.env)
}

#-------------------------------------------------------------------------------
#' Optimize a proliferation model using maximum likelihood
#'
#' @param x numeric vector of raw log10 intensities (one value per event).
#' @param starts named list of starting parameter values.
#' @param upper named numeric vector of upper bounds for each parameter.
#' @param lower named numeric vector of lower bounds for each parameter.
#' @param fixed named numeric vector of fixed parameters excluded from optimisation.
#' @param peak.stats data frame of initial peak estimates from \code{find_initial_peaks}.
#' @param plot logical; if TRUE produce diagnostic plots of the fit (default FALSE).
#' @param opt.env optional environment for logging parameter traces (default NULL).
#' @param verbose logical; if TRUE print parameter values at each iteration (default FALSE).
#' @param log logical; if TRUE use log-scale density in the likelihood (default TRUE).
#' @returns named numeric vector of optimised parameters (free + fixed).
opt_prolif_model_mle <- function(x, starts, upper, lower, fixed, peak.stats, plot=F, opt.env=NULL, verbose=F, log=T) {

  res <- optim(
    par = starts,
    fn = prolif_model_wrapper,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    log=log,
    type="neg_log_likelihood",
    method="L-BFGS-B",
    verbose=verbose
  )

  final.par <- c(res$par, fixed)

  if (plot) {

    if (!is.null(opt.env)) {
      opt_plot_params(x, opt.env)
    }

    hd <- hist(x, breaks=250, plot=F)
    plot(hd$mids,
         hd$density,
         xlab="Log10(intensity)",
         ylab="Density",
         pch=20,
         col="grey", bty="n",
         main="Estimated fit vs optimal fit")

    mod.d <- prolif_model_wrapper(starts, nrow(peak.stats), x = hd$mids, type="density")
    lines(hd$mids,
          exp(mod.d),
          col="blue")

    mod.d <- prolif_model_wrapper(final.par, nrow(peak.stats), x = hd$mids, type="density")
    lines(hd$mids,
          exp(mod.d),
          col="red")

    legend("topleft",
           legend=c("Initial fit", "Optimized fit", "Data"),
           fill=c("blue", "red", "grey"),
           bty="n")
  }

  return(final.par)
}

#-------------------------------------------------------------------------------
#' Optimize a proliferation model using least squares
#'
#' @param x numeric vector of histogram midpoints (log10 scale).
#' @param y numeric vector of smoothed counts matching x in length.
#' @param starts named list of starting parameter values.
#' @param upper named numeric vector of upper bounds for each parameter.
#' @param lower named numeric vector of lower bounds for each parameter.
#' @param fixed named numeric vector of fixed parameters excluded from optimisation.
#' @param peak.stats data frame of initial peak estimates from \code{find_initial_peaks}.
#' @param plot logical; if TRUE produce diagnostic plots of the fit (default FALSE).
#' @param opt.env optional environment for logging parameter traces (default NULL).
#' @returns named numeric vector of optimised parameters (free + fixed).
opt_prolif_model_ls <- function(x, y, starts, upper, lower, fixed, peak.stats, plot=F, opt.env=NULL) {

  res <- minpack.lm::nls.lm(
    par = starts,
    fn = prolif_resid,
    y = y,
    x = x,
    fixed=fixed,
    n.peaks = nrow(peak.stats),
    upper=upper,
    lower=lower,
    opt.env=opt.env,
    control = minpack.lm::nls.lm.control(
      nprint = F,
      maxiter = 1024,
      factor = 0.01,
      maxfev=1000000,
      ptol=.Machine$double.xmin,
      gtol= 0,
      ftol=.Machine$double.xmin)
  )

  final.par <- c(res$par, fixed)

  if (!is.null(opt.env)) {
    opt.env[["score"]]     <- log2(res$rsstrace)
  }

  if (plot) {

    if (!is.null(opt.env)) {
      opt_plot_params(x, opt.env)
    }

    plot(x, y,
         main="Estimated fit vs optimal fit",
         xlab="Bin average log10(CTV)",
         ylab="Count",
         pch=20,
         col="grey", bty="n")

    y.pred <- sapply(x, function(x){prolif_model(means=peak.stats$est_mean,
                                                  sd=peak.stats$est_peak_sd,
                                                  summits=peak.stats$est_summit,
                                                  x)})

    lines(y.pred ~ x, lwd=2, col="blue")

    y.pred.optim <- sapply(x, function(x){prolif_model_wrapper(final.par, nrow(peak.stats), x=x, fixed=fixed)})

    lines(y.pred.optim ~ x, lwd=2, col="red")
    legend("topleft",
           legend=c("Initial fit", "Optimized fit", "Smoothed data"),
           fill=c("blue", "red", "grey"),
           bty="n")
  }

  return(final.par)
}


#-------------------------------------------------------------------------------
#' Fit a proliferation model using least squares or maximum likelihood
#'
#' Fit peaks on a CTV or CFSE or other tracking dye trace using a binned
#' or maximum likelihood approach. The trace should be raw per event data from
#' the FCS file on the CTV+ population.
#'
#' @param trace raw FACS intensities
#' @param peak.0.lower.bound the value on log10 scale where the first valley is
#' @param peak.thresh.enrich fold change over valley to call peak in initial
#' estimation (default 1)
#' @param peak.thresh.summit minimum height of a peak in percentage of total
#' heights (default 0.05)
#' @param peak.max the maximum number of peaks to search for (default 12)
#' @param bins number of bins for fitting (default 250)
#' @param smoothing.window number of values up and downstream for the NN
#' smoother (default 2)
#' @param window.scaling.factors two-element vector scaling the search window
#' left and right of each expected peak position during initial estimation
#' (default c(0.25, 0.25)). See details.
#' @param plot logical; if TRUE produce all plots. Overrides plot.optim and
#' plot.final (default TRUE).
#' @param plot.optim logical; if TRUE produce optimisation diagnostic plots
#' (default TRUE).
#' @param plot.final logical; if TRUE produce final per-peak fit plot (default TRUE).
#' @param plot.main character string title for the plot (default "Proliferation model").
#' @param opt.peak.pos.dev numeric; maximum allowed deviation of each peak mean
#' from its initial estimate during optimisation. Defaults to half the estimated
#' inter-peak distance (NULL). Set to Inf to remove constraint. See details.
#' @param opt.trim.left.tail logical; if TRUE remove data points left of the
#' leftmost initial peak estimate before optimisation (default FALSE). See details.
#' @param opt.trim.right.tail logical; if TRUE remove data points right of the
#' rightmost initial peak estimate before optimisation (default FALSE). See details.
#' @param mode character string; \code{"LS"} for Levenberg-Marquardt non-linear
#' least squares (default) or \code{"MLE"} for maximum likelihood estimation.
#' See details.
#' @param full.out logical; if TRUE return a list containing peak statistics,
#' histogram object, x values, model predictions, and final parameters. If FALSE
#' return only the peak statistics data frame (default FALSE).
#' @param peak.x.model logical; if TRUE fit an additional unconstrained peak
#' (peak X) to capture a dye-negative cell population (default FALSE). See details.
#' @param peak.x.upper.bound numeric; hard upper limit (log10 scale) for the
#' position of peak X. Set to NULL for no limit (default NULL).
#' @param peak.x.position numeric; manually specify the starting position of
#' peak X (log10 scale). If NULL, position is estimated automatically (default NULL).
#' @param peak.x.fixed logical; if TRUE fix peak X at \code{peak.x.position}
#' and do not optimise its mean (default FALSE).
#' @param peak.x.thresh.summit numeric; minimum height of peak X relative to the
#' trace mode for autodetection. Values between 0 and 1 (default 0.05).
#' @param peak.x.thresh.enrich numeric; minimum fold enrichment of peak X over
#' its flanking valleys for autodetection (default 1.2).
#' @param peak.x.sd numeric; upper bound for the SD of peak X during
#' optimisation. If NULL, the upper bound defaults to the same as the shared
#' \code{peak.sd} upper bound (default NULL).
#' @param verbose logical; if TRUE print optimisation parameter values at each
#' iteration (default FALSE).
#' @param log logical; if TRUE use log-scale density in MLE mode (default TRUE).
#' @param ... additional arguments passed to \code{find_initial_peaks}.
#'
#' @details
#' This fits a proliferation model on a trace of raw FACS intensities.
#' The fitting happens in two stages:
#' 1. Initial estimation of number of peaks and their positions
#' 2. Optimization of the model to find final values
#'
#' Initial estimation is done based on thresholds. First the 0 peak
#' is estimated by finding the mode of the data > peak.0.lower.bound.
#'
#' The next peak is then found by finding the mode, at half of the intenstiy
#' of the mode of peak 0. The range where the mode is looked for can be
#' controlled by tweaking `window.scaling.factors`.
#'
#' Peaks are only called if they exceed `peak.thresh.enrich`, the relative
#' enrichment between the valley's either side of the peak and the peak summit
#' and `peak.thresh.summit`, the height of the peak with respect to the mode of
#' the whole trace. So this takes a value between 0-1 describing a percentage of
#' the mode the peak must have. To fit all possible peaks, set both of these to 0
#'
#' Peaks must also be adjacent. So if peak 2 does not pass `peak.thresh.enrich`
#' but peak 3 does, only 2 peaks are fit. This is done to avoid fitting many
#' peaks in the marginal count range.
#'
#' After estimating starting values, a Gaussian mixture distribution is fit to
#' the trace (a sum of the PDF of individuals Gaussian). The starting values
#' for the peak standard deviations are estimated from the trace in the bound
#' of peak0. These values are input into the Levenberg-Marquardt algo to
#' optimize with respect to the residual sum of squares between the model
#' and the smoothed trace in mode 'LS' or using MLE and R optim() when in mode
#' 'MLE'.
#'
#' # `opt.peak.pos.dev`
#' opt.peak.pos.dev controls the range peaks are allowed to vary in position
#' during optimization. Defaults to +- half the estimated peak distance (NULL).
#' To ignore this, set to Inf. This setting avoids the optimizer putting peaks
#' in the left tail that sometimes might be present or separating peaks into
#' nonsensical distances. To contol the allowable range of initial peak
#' estimates see `window.scaling.factors`
#'
#' # `opt.trim.left.tail` / `opt.trim.right.tail`
#' Removes values (MLE) or sets count to zero (LS) of values that fall outside
#' 1sd of the initial model space. This essentially functions as an auto gate
#' after identifying initial peak positions.
#'
#' # `mode`
#' The package offers two optimization schemes:
#'
#' Non linear least squares (LS) with mode="LS". Here the data is first binned
#' and smoothed, then parameters are optimized over the smoothed trace using non
#' linear least squares. This method is quick and robust, but relies on binning
#' and some other data processing so is technically less accurate.
#'
#' Maximum likelihood estimation (MLE) based with mode="MLE". This mode is a
#' little more proper in that it does not rely on binning or other data tricks
#' to optimize, and takes the full data set to optimize on directly. Initial
#' parameter estimation is still done on a binned smoothed version of the data,
#' but optimization is not.
#'
#' In practice I have found very little difference between these two when the
#' setup is performed correctly.
#'
#' # `window.scaling.factors`
#' window.scaling.factors control the size of the window left and right around the next peak position
#' to find the mode of the next peak. The smaller this value is, the closer to the half intensity the
#' peak estimates will be.
#'
#' I.e. if the estimated distance between peaks is 0.5 and the scaling
#' factors are c(0.25, 0.5) and the current peak is 1.5 the next peak mode will be estimates as:
#' lower <- 0.5 * 0.25
#' upper <- 0.5 * 0.5
#'
#' # Where should the next peak be based on half it's intensity
#' est.peak.pos <- log10(10^1.5/2)
#'
#' Then find the max value in the window est.peak.pos-lower, est.peak.pos+upper
#'
#' @returns A data frame with peak statistics or list of objects if full.out=T
#'
#' @examples
#'
#' # Simulate proliferation data
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c(y/4, y/2, y)
#'
#' # Fit peaks
#' peaks <- fit_peaks(y, peak.0.lower.bound=9.8)
#'
#' # Simulate proliferation data #2
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c((y/8)[1:500], y/4, y/2, y[1:500])
#'
#' peaks <- fit_peaks(y, peak.0.lower.bound=9.9, mode="MLE")
#'
#' @export
fit_peaks  <- function(trace,
                       peak.0.lower.bound,
                       bins=250,
                       smoothing.window=2,
                       plot.optim=T,
                       plot.final=T,
                       plot=T,
                       plot.main = "Proliferation model",
                       opt.peak.pos.dev=NULL,
                       window.scaling.factors=c(0.25,0.25),
                       opt.trim.left.tail=F,
                       opt.trim.right.tail=F,
                       mode="LS",
                       full.out=F,
                       peak.x.model=F,
                       peak.x.upper.bound=NULL,
                       peak.x.position=NULL,
                       peak.x.fixed=F,
                       peak.x.thresh.summit=0.05,
                       peak.x.thresh.enrich=1.2,
                       peak.x.sd=NULL,
                       verbose=F,
                       log=T,
                       ...) {

  # Plotting
  if (!plot) {
    plot.optim <- F
    plot.final <- F
  }
  if (plot.optim) {par(mfrow=c(2,3), mar=c(5,5,5,1))}

  # Construct histogram
  cur.hist <- hist(log10(trace), breaks=bins, plot=F)

  # Nearest neighbor smoother
  y.smth     <- nn_smoother(cur.hist$counts, window=smoothing.window)
  y.raw      <- log10(trace)

  # Find the initial peak estimates
  peak.stats <- find_initial_peaks(cur.hist$mids, y.smth, peak.0.lower.bound, plot = plot.optim, window.scaling.factors=window.scaling.factors,...)

  #-----------------------------------------------------------------------------
  # Roughly estimate peak X and update stats table
  if (peak.x.model) {
    if (!peak.x.fixed) {

      if (!is.null(peak.x.position)) {
        peak.x.mode <- peak.x.position
      } else {
        peak.x.mode <- find_peak_x_approx_mode(y.raw,
                                               peak.x.thresh.summit=peak.x.thresh.summit,
                                               peak.x.thresh.enrich=peak.x.thresh.enrich,
                                               peak.x.upper.bound=peak.x.upper.bound)
      }

      if (!is.null(peak.x.mode)) {

          est.smt       <- y.smth[nearest_index(cur.hist$mids, peak.x.mode)]
          est.peak.dist <- peak.stats[1,1] - peak.stats[2,1]
          est.enrich    <- find_enrichment(cur.hist$mids, y.smth, peak.x.mode, est.peak.dist/2)

          peak.stats    <- peak.stats[peak.stats$est_mean > peak.x.mode,]

          peak.stats    <- rbind(peak.stats,
                                 c(peak.x.mode, est.smt, est.smt, nrow(peak.stats)+1, nrow(peak.stats), 0))

          peak.stats$est_summit_percentage <- peak.stats$est_summit / sum(peak.stats$est_summit)

      } else {
        peak.x.mode = NULL
        peak.x.model = FALSE
        warning("No valid mode for peak.x found, skipping fitting. Adjust parameters if it should be there.")
      }

    } else {
      peak.x.mode = peak.x.position
      peak.x.model = TRUE
    }
  } else {
    peak.x.mode = NULL
    peak.x.model = FALSE
  }

  # Plot
  if (plot.optim) {
   opt_plot_estimates(cur.hist, y.smth, peak.stats, peak.0.lower.bound, peak.x.model, peak.x.mode, peak.x.upper.bound)
  }

  #-----------------------------------------------------------------------------
  # Estimate starting parameters
  #-----------------------------------------------------------------------------
  x.mids        <- cur.hist$mids
  est.gen0.mean <- peak.stats[1, "est_mean"]
  est.gen1.mean <- peak.stats[2, "est_mean"]

  if (is.na(est.gen1.mean)) {
    est.gen1.mean <- log10((10^est.gen0.mean)/2)
  }

  est.peak.dist <- est.gen0.mean - est.gen1.mean

  #-----------------------------------------------------------------------------
  if (opt.trim.left.tail) {
    lower.lim                  <- peak.stats[nrow(peak.stats),"est_mean"] -  (est.peak.dist/2)
    y.smth[x.mids < lower.lim] <- 0
    y.raw                      <- y.raw[y.raw > lower.lim]
  }

  if (opt.trim.right.tail) {
    upper.lim                  <- peak.stats[1,"est_mean"] +  (est.peak.dist/2)
    y.smth[x.mids > upper.lim] <- 0
    y.raw                      <- y.raw[y.raw < upper.lim]
  }

  #-----------------------------------------------------------------------------
  # Roughly estimate the SD of generation zero
  peak0.mean  <- peak.stats[1,"est_mean"]
  h.dist      <- est.peak.dist / 2
  yy          <- y.smth[(x.mids > est.gen0.mean-h.dist) & (x.mids < est.gen0.mean + h.dist)]
  xx          <- x.mids[(x.mids > peak0.mean-h.dist) & (x.mids < peak0.mean+h.dist)]
  est.peak.sd <- sd(xx)
  peak.stats[,"est_peak_sd"] <- est.peak.sd

  #-----------------------------------------------------------------------------
  # Starting parameter list
  starts          <- c(as.list(peak.stats$est_mean), as.list(peak.stats$est_summit))
  mean.names      <- paste0("gen", 0:(nrow(peak.stats)-1), ".mean")
  summit.names    <- paste0("gen", 0:(nrow(peak.stats)-1), ".summit")
  names(starts)   <- c(mean.names, summit.names)
  starts$peak.sd  <- est.peak.sd

  if (peak.x.model) {
    sd.names <- c(rep("peak.sd", nrow(peak.stats)-1), "genX.sd")
  } else {
    sd.names <- rep("peak.sd", nrow(peak.stats))
  }

  fixed  <- NULL
  starts <- starts[!names(starts) %in% names(fixed)]

  #-----------------------------------------------------------------------------
  # Determine the limits of optimization of peak positions
  upper        <- rep(Inf, length(starts))
  names(upper) <- names(starts)
  lower        <- rep(0, length(starts))
  names(lower) <- names(starts)

  if (is.null(opt.peak.pos.dev)) {
    opt.peak.pos.dev <- est.peak.dist/2
  }

  par.means        <- grep("mean", names(starts))
  upper[par.means] <- as.numeric(starts[par.means]) + (opt.peak.pos.dev * window.scaling.factors[2])
  lower[par.means] <- as.numeric(starts[par.means]) - (opt.peak.pos.dev * window.scaling.factors[1])

  lower["peak.sd"] <-1e-16

  # Create environment to store iterations
  opt.env <- opt_new_env(nrow(peak.stats))

  # Add the constraints and start for genX
  if (peak.x.model) {
    starts[["genX.sd"]] <- starts$peak.sd*2
    lower[["genX.sd"]]  <- lower["peak.sd"]
    upper[["genX.sd"]]  <- if (!is.null(peak.x.sd)) peak.x.sd else upper["peak.sd"]

    if (!is.null(peak.x.upper.bound)) {
      upper[[par.means[length(par.means)]]] <- peak.x.upper.bound
    } else {
      upper[[par.means[length(par.means)]]] <- upper[[par.means[length(par.means)]]]+starts$peak.sd*2
    }

    if (peak.x.fixed) {
      if (!is.null(peak.x.position)) {
        if (verbose) {cat("[INFO] Fixing peak x position")}
        starts[[par.means[length(par.means)]]] <- peak.x.position
        upper[[par.means[length(par.means)]]] <- peak.x.position
        lower[[par.means[length(par.means)]]] <- peak.x.position
      } else {
        stop("peak.x.fixed=T but no position provided")
      }
    }
  }

  #-----------------------------------------------------------------------------

  if (verbose) {
    cat("[INFO] starting params: ", "\n")
    print(starts)
    cat("[INFO] upper:  ", "\n")
    print(upper)
    cat("[INFO] lower:  ", "\n")
    print(lower)
    cat("[INFO] number of bins:  ", length(x.mids), "\n")
  }

  #-----------------------------------------------------------------------------
  # Run optimization
  if (mode == "LS") {
    final.par <- opt_prolif_model_ls(x.mids,
                                     y.smth,
                                     starts=starts,
                                     upper=upper,
                                     lower=lower,
                                     fixed=fixed,
                                     peak.stats=peak.stats,
                                     plot=plot.optim,
                                     opt.env=opt.env)

    y.pred.optim <- prolif_model_wrapper(final.par,
                                         nrow(peak.stats),
                                         x=x.mids)

    total.events <- sum(y.smth)

    peak.stats$opt_summit        <- as.numeric(final.par[summit.names])

    for (i in 1:nrow(peak.stats)) {
      int.res <- integrate(prolif_single_peak,
                           lower=min(y.raw),
                           upper=max(y.raw),
                           mean=as.numeric(final.par[mean.names[i]]),
                           summit=as.numeric(final.par[summit.names[i]]),
                           sd=as.numeric(final.par[sd.names[i]]))

      peak.stats[i,"peak_area"]       <- int.res$value
      peak.stats[i,"peak_area_error"] <- int.res$abs.error
    }

  } else if (mode == "MLE") {
    starts[summit.names] <- as.numeric(starts[summit.names]) / sum(as.numeric(starts[summit.names]))
    upper[summit.names]  <- 1

    final.par <- opt_prolif_model_mle(y.raw,
                                      starts=starts,
                                      upper=upper,
                                      lower=lower,
                                      fixed=fixed,
                                      peak.stats=peak.stats,
                                      plot=plot.optim,
                                      opt.env=opt.env,
                                      verbose=verbose,
                                      log=log)

    total.events <- length(y.raw)
    binwidth     <- (x.mids[2] - x.mids[1])

    y.pred.optim <- exp(prolif_model_wrapper(final.par,
                                             nrow(peak.stats),
                                             x=x.mids,
                                             type="density"))

    y.pred.optim <- (y.pred.optim*binwidth) * total.events

    y.pred.summit <- exp(prolif_model_wrapper(final.par,
                                              nrow(peak.stats),
                                              x=as.numeric(final.par[mean.names]),
                                              type="density"))
    y.pred.summit <- (y.pred.summit*binwidth) * total.events

    peak.stats$opt_summit        <- y.pred.summit
    peak.stats$peak_area         <- final.par[summit.names]
    peak.stats$peak_area_error   <- 0
  } else {
    stop(simpleError("No valid mode provided, must be LS, MIXDIST or MLE."))
  }

  peak.stats$opt_mean          <- as.numeric(final.par[mean.names])
  peak.stats$opt_peak_sd       <- as.numeric(final.par[sd.names])

  if (plot.optim) {
    mtext(plot.main, side = 3, line = -2, outer = TRUE)
    par(mfrow=c(1,1))
  }

  rm("opt.env")

  #-----------------------------------------------------------------------------
  peak.stats$peak_area_prop   <- (peak.stats$peak_area / sum(peak.stats$peak_area)) *100
  peak.stats$peak_summit_prop <- (peak.stats$opt_summit / sum(peak.stats$opt_summit)) *100
  peak.stats$peak_events      <- total.events * (peak.stats$peak_area_prop/100)

  tmp <- vector()
  for(i in 1:nrow(peak.stats)) {
    tmp[i] <- (peak.stats[i, "peak_events"] / 2^(i-1))
  }
  peak.stats$peak_ancestors <- tmp

  if (peak.x.model) {
    peak.stats[nrow(peak.stats), "generation"] <- "x"
  }

  #-----------------------------------------------------------------------------
  if (plot.final) {
    opt_plot_final_pp(x.mids, cur.hist$counts, peak.stats, main=plot.main,
                      peak.x.model=peak.x.model,
                      peak.0.lower.bound=peak.0.lower.bound,
                      peak.x.upper.bound=peak.x.upper.bound)
  }

  if (full.out) {
    return(list(peak.stats=peak.stats,
                hist=cur.hist,
                x=x.mids,
                y=y.pred.optim,
                par=final.par))
  } else {
    return(peak.stats)
  }
}

#-------------------------------------------------------------------------------
#' Get proliferation statistics
#'
#' Takes output from fit_peaks and calculates proliferation statistics.
#' Proliferation statistics according to the FCS express definitions
#' https://fcsexpressdownloads.s3.amazonaws.com/manual/manual_WIN_RUO/index.html?proliferation_statistics.htm
#'
#' Can be used manually by provided a data frame with two columns
#'
#' - generation
#'
#' - peak_events
#'
#' Where each row represent the peak number and number of events in a peak.
#'
#' @param peak.stats output from fit_peaks (data frame)
#' @returns A dataframe with proliferation statistics
#' @examples
#'
#' # Simulate proliferation data
#' y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
#' y <- c(y/4, y/2, y)
#'
#' # Fit peaks
#' peaks <- fit_peaks(y, peak.0.lower.bound=9.8)
#' get_prolif_stats(peaks)
#'
#' # Manual table
#' peaks <- data.frame(generation=c(0,1,2), peak_events=c(1000, 1000, 1000))
#' get_prolif_stats(peaks)
#' @export
get_prolif_stats <- function(peak.stats) {

  total.events <- sum(peak.stats[,"peak_events"])
  gen.x.events <- peak.stats[peak.stats$generation=="x","peak_events"]

  if (length(gen.x.events) == 0) {
    gen.x.events <- NA
  }

  peak.stats <- peak.stats[peak.stats$generation != "x",]
  peak.stats$generation <- as.numeric(peak.stats$generation)

  founding.pop <- sum(peak.stats$peak_events / 2^peak.stats$generation)
  divided.pop  <- sum(peak.stats$peak_events[-1] / 2^peak.stats$generation[-1])

  prol.index   <- sum(peak.stats$peak_events) / founding.pop
  div.index    <- sum(peak.stats$peak_events[-1]) / divided.pop
  perc.div     <- (divided.pop / founding.pop)*100

  vec <- c(founding.pop, divided.pop, prol.index, div.index, perc.div, nrow(peak.stats)-1, total.events, gen.x.events)
  names(vec) <- c("founding_population", "divided_population", "proliferation_index", "division_index", "percent_divided", "num_generations", "total_events", "events_peak_x")

  return(vec)
}
