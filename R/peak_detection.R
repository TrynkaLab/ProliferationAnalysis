#-------------------------------------------------------------------------------
#' Find the next peak in a trace
#'
#' @description
#' Find the mode of the next peak using the previous peak and the average peak
#' distance. The window can be asymmetrically scaled using the scaling factors.
#' This is done as the smaller the peaks get, the closer together they tend to be.
find_next_peak <- function(x, y, prev.peak, peak.dist, window.scaling.factors=c(0.25, 0.25)) {

  est.curpeak     <- log10((10^prev.peak)/2)

  low  <- est.curpeak - (peak.dist * window.scaling.factors[1])
  high <- est.curpeak + (peak.dist * window.scaling.factors[2])

  curpeak.mean    <- find_local_max(x, y, low, high)
  curpeak.summit  <- y[nearest_index(x, curpeak.mean)]

  fold.change  <- find_enrichment(x, y, curpeak.mean, peak.dist/2)

  return(c(curpeak.mean, curpeak.summit, fold.change))
}

#-------------------------------------------------------------------------------
#' Find the relative enrichment over a valley
find_enrichment <- function(x, y, curpeak.mean, est.window) {

  curpeak.summit  <- y[nearest_index(x, curpeak.mean)]

  a <- nearest_index(x, curpeak.mean + est.window)
  b <- nearest_index(x, curpeak.mean - est.window)
  valley.count <- mean(y[c(a, b)])
  fold.change  <- curpeak.summit/valley.count

  return(fold.change)
}

#-------------------------------------------------------------------------------
#' Estimate initial parameters for proliferation model
#'
#' @param x x values (histogram midpoints)
#' @param y a smoothed or raw count trace
#' @param peak.0.lower.bound lower bound for finding gen0
#' @param peak.thresh.enrich fold change over valley to call peak in initial estimation (default 1)
#' @param peak.thresh.summit minimum height of a peak in percentage of total heights (default 0.05)
#' @param peak.max the maximum number of peaks to search for (default 12)
#' @param plot should plot be generated
#' @param window.scaling.factors Scaling factors for peak differences. See details (default c(0.25, 0.25))
#'
#' @details
#' `window.scaling.factors`
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
#' @returns data frame with initial peak estimates
find_initial_peaks <- function(x, y, peak.0.lower.bound, peak.thresh.enrich=1, peak.thresh.summit=0.05, peak.max=12, plot=F, window.scaling.factors=c(0.25,0.25)) {
  peak0.mean       <- find_local_max(x, y, peak.0.lower.bound, max(x))

  if (is.na(peak0.mean)) {
    cat("[WARN] No values found in peak.0.lower.bound. Returning NA.\n")
    return(NULL)
  }

  peak0.summit     <- y[nearest_index(x, peak0.mean)]

  est.peak1     <- log10(10^peak0.mean/2)
  est.peak.dist <- peak0.mean - est.peak1
  peak0.enrichment <- find_enrichment(x, y, peak0.mean, est.peak.dist/2)

  prev.peak.mean <- peak0.mean
  peak.stats     <- data.frame()
  lim            <- min(x)

  for (i in 1:peak.max) {

    if (prev.peak.mean-est.peak.dist < lim) {
      next;
    }
    res            <- find_next_peak(x, y, prev.peak.mean, est.peak.dist, window.scaling.factors=window.scaling.factors)
    peak.stats     <- rbind(peak.stats, res)
    prev.peak.mean <- res[1]

  }

  peak.stats[is.na(peak.stats)] <- 0

  peak.stats            <- rbind(c(peak0.mean,
                                   peak0.summit,
                                   peak0.enrichment),
                                 peak.stats)
  colnames(peak.stats)  <- c("est_mean", "est_summit", "est_enrichment")
  peak.stats$index      <- 1:(nrow(peak.stats))
  peak.stats$generation <- 0:(nrow(peak.stats)-1)

  peak.stats$est_summit_percentage <- peak.stats$est_summit/sum(peak.stats$est_summit)

  j <- 1
  for (i in 2:nrow(peak.stats)) {

    if (peak.stats[i, "est_summit_percentage"] > peak.thresh.summit) {
      j <- j+1
      next
    } else {
      if (peak.stats[min(i+1, nrow(peak.stats)), "est_summit_percentage"] > peak.thresh.summit) {
        j <- j+1
        next
      } else {
        break
      }
    }
  }
  peak.stats <- peak.stats[1:j,]

  peak.stats <- peak.stats[(peak.stats$est_enrichment > peak.thresh.enrich) | peak.stats$generation == 0,, drop=F]

  n.peaks  <- 0
  for (i in 1:(nrow(peak.stats)-1)) {
    cur.idx <- peak.stats[i, "index"]

    if (cur.idx==1) {
      n.peaks  <- 1
    }

    if (i+1 > nrow(peak.stats)) {
      break
    }

    if (cur.idx +1 == peak.stats[i+1, "index"]){
      n.peaks  <- n.peaks +1
    }  else {
      break
    }
  }

  if (n.peaks==0) {
    n.peaks <- 1
  }

  peak.stats <- peak.stats[1:n.peaks,, drop=F]

  return(peak.stats)
}

#-------------------------------------------------------------------------------
#' Find the mode of the first left sided peak in the data, given some constraints.
#'
#' @param trace raw FACS intensities (log10 scale)
#' @param peak.x.thresh.summit minimum relative height of peak x
#' @param peak.x.thresh.enrich minimum enrichment of peak x over valley
#' @param peak.x.upper.bound hard upper limit for peak x position
#'
#' @export
find_peak_x_approx_mode <- function(trace, peak.x.thresh.summit, peak.x.thresh.enrich, peak.x.upper.bound=NULL, ...) {

  if (is.null(peak.x.upper.bound)) {
    peak.x.upper.bound <- max(trace)
  }

  dens <- density(trace, ...)

  d     <- dens$y / max(dens$y)

  dd    <- d[d > peak.x.thresh.summit]

  if (length(dd) < 3) {
    return(NULL)
  }

  xx    <- dens$x[d >peak.x.thresh.summit]

  max.d <- 0
  for (i in 1:length(dd)) {
    if (dd[i] > max.d) {
      max.d <- dd[i]
    }
    if ((max.d/dd[i]) >= peak.x.thresh.enrich) {
      break
    }
  }

  max.x <- xx[which(dd == max.d)]

  if (max.x > peak.x.upper.bound) {
    msg <- "No peak found at specified thresholds. Returning the mode witin hard limit set by peak.x.upper.bound"
    warning(simpleWarning(msg))

    return(NULL)
    # TODO: Check if this can be overidden
    tmp <- dd[xx < peak.x.upper.bound]

    if (length(tmp) > 2) {
      return(xx[which(dd ==  max(tmp))])
    } else {
      return(NULL)
    }

  } else {
    return(max.x)
  }
}
