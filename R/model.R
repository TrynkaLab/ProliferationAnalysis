#-------------------------------------------------------------------------------
#' Scaled distribution of a single peak in prolif model
#'
#' @param x x values (log scaled intensity values)
#' @param mean the mean of the peak
#' @param sd the sd of the peak
#' @param summit the peak height (estimated relative to the other peaks)
#' @returns scaled density for peak height
prolif_single_peak <- function(x, mean, sd, summit) {
  cur.density <- dnorm(x, mean=mean, sd=sd)
  cur.density <- cur.density / dnorm(mean, mean=mean, sd=sd)
  return(cur.density*summit)
}

#-------------------------------------------------------------------------------
#' Mixture distribution of Gaussian
#'
#' @description
#' Gaussian prolif model, it is a mixture distribution scaled by relative
#' peak heights.
#'
#' @param means numeric vector of peak means (log10 scale).
#' @param sd numeric vector of peak standard deviations, one per peak.
#' @param summits numeric vector of peak heights (counts or relative weights).
#' @param x numeric vector of x values at which to evaluate the model.
#' @returns numeric vector of predicted y values (sum of scaled Gaussian densities).
prolif_model <- function(means, sd, summits, x) {

  y.pred <- 0

  for (i in 1:length(means)) {
    cur.pred <- prolif_single_peak(x, means[i], sd[i], summits[i])
    y.pred   <- y.pred + cur.pred
  }

  return(y.pred)
}

#-------------------------------------------------------------------------------
#' Density of mixture distribution of Gaussian
#'
#' @param means numeric vector of peak means (log10 scale).
#' @param sd numeric vector of peak standard deviations, one per peak.
#' @param summits numeric vector of mixture weights (will be normalised to sum to 1).
#' @param x numeric vector of x values at which to evaluate the density.
#' @param verbose logical; if TRUE print parameter values at each call (default FALSE).
#' @param log logical; if TRUE return log-density using the log-sum-exp trick for
#'   numerical stability (default FALSE).
#' @param opt.env optional environment created by \code{opt_new_env} for logging
#'   parameter traces during optimisation (default NULL).
#' @returns numeric vector of density (or log-density) values at each x.
prolif_model_density <- function(means, sd, summits, x, verbose=F, log=F, opt.env=NULL) {
  density <- 0
  summits <- summits/sum(summits)

  if (verbose) {
    cat("[INFO] means: ",   means, "\n")
    cat("[INFO] sd: ",      sd, "\n")
    cat("[INFO] summits: ", summits, "\n")
    cat("[INFO] --------------------------------------------- \n")
  }

  if (log){
    # https://stats.stackexchange.com/questions/105602/example-of-how-the-log-sum-exp-trick-works-in-naive-bayes
    densities <- matrix(nrow=length(means), ncol=length(x))
    for (i in 1:length(means)) {
      cur.d         <- dnorm(x, mean=means[i], sd=sd[i], log=T)
      densities[i,] <- log(summits[i]) + cur.d
    }

    density <- apply(densities, 2, function(cur.dens) {
      a        <- max(cur.dens)
      cur.dens <- exp(cur.dens - a)
      return(a + log(sum(cur.dens)))
    })

  } else {
    density <- 0
    for (i in 1:length(means)) {
      cur.d   <- dnorm(x, mean=means[i], sd=sd[i])
      density  <- density + (summits[i] * cur.d)
    }

    # Hack to deal with numeric precision limit on very small densities
    if (sum(density==0) >0) {
      msg <- "[DO NOT IGNORE] Zero densities detected. Consider using log=T. This is due to a numeric precision limit in dnorm().
    As a hack setting zeroes to 5e-324. This does invalidate the PDF!!"
      warning(simpleWarning(msg))
      density[density==0] <- 5e-324
    }
  }

  return(density)
}

#-------------------------------------------------------------------------------
#' Negative log likelihood of mixture distribution of Gaussian
#'
#' @param means numeric vector of peak means (log10 scale).
#' @param sd numeric vector of peak standard deviations, one per peak.
#' @param summits numeric vector of mixture weights (will be normalised to sum to 1).
#' @param x numeric vector of observed data values (log10 intensities).
#' @param verbose logical; if TRUE print the NLL at each call (default FALSE).
#' @param log logical; if TRUE compute density on the log scale for numerical
#'   stability (default TRUE).
#' @param invert logical; if TRUE return the positive log likelihood instead of
#'   the negative (default FALSE).
#' @param opt.env optional environment created by \code{opt_new_env} for logging
#'   parameter traces during optimisation (default NULL).
#' @returns scalar negative log likelihood value.
prolif_model_nll <- function(means, sd, summits, x, verbose=F, log=T, invert=F, opt.env=NULL) {
  density <- prolif_model_density(means, sd, summits, x, verbose=verbose, log=log, opt.env=opt.env)

  if (log) {
    nll <- -sum(density)
  } else {
    nll <- -sum(log(density))
  }

  if (invert) {
    nll <- -nll
  }

  if (verbose) {
    cat("[INFO] nll: ", nll, "\n")
  }
  return(nll)
}

#-------------------------------------------------------------------------------
#' Wrapper that returns a gaussian proliferation model.
#'
#' @description
#' Takes a named parameter vector, unpacks it into means, summits, and SDs, and
#' dispatches to the requested model type. Used as the objective function passed
#' to \code{optim} or \code{nls.lm}.
#'
#' @param par named numeric vector of model parameters. Expected names follow the
#'   convention \code{gen0.mean}, \code{gen0.summit}, ..., \code{peak.sd} and
#'   optionally \code{genX.sd}.
#' @param n.peaks integer number of peaks in the model.
#' @param x numeric vector of x values (histogram midpoints or raw log10 intensities).
#' @param fixed named numeric vector of fixed parameters to append to \code{par}
#'   before unpacking (default NULL).
#' @param type character string selecting the return value: \code{"prolif_model"}
#'   (default, returns scaled counts), \code{"density"} (returns log-density),
#'   or \code{"neg_log_likelihood"} (returns NLL scalar).
#' @param verbose logical; if TRUE print parameter values at each call (default FALSE).
#' @param opt.env optional environment created by \code{opt_new_env} for logging
#'   parameter traces during optimisation (default NULL).
#' @param invert logical; passed to \code{prolif_model_nll} when
#'   \code{type="neg_log_likelihood"} (default FALSE).
#' @param names character vector of parameter names; overridden by \code{names(par)}
#'   if present (default NULL).
#' @param log logical; if TRUE use log-scale density computation (default TRUE).
#' @returns numeric vector (model predictions or density) or scalar NLL,
#'   depending on \code{type}.
prolif_model_wrapper <- function(par, n.peaks, x, fixed=NULL, type="prolif_model", verbose=F, opt.env=NULL, invert=F, names=NULL, log=T) {

  if (!is.null(names(par))) {
    names      <- names(par)
  }
  par        <- as.numeric(par)
  names(par) <- names

  mean.names   <- paste0("gen", 0:(n.peaks-1), ".mean")
  summit.names <- paste0("gen", 0:(n.peaks-1), ".summit")

  if (!is.null(fixed)) {
    par <- c(par, fixed)
  }

  mean.vec   <- as.numeric(par[mean.names])
  summit.vec <- as.numeric(par[summit.names])
  sd.vec     <- as.numeric(par["peak.sd"])

  if (length(sd.vec)==1) {
    sd.vec <- rep(sd.vec, length(mean.vec))
  }

  if ("genX.sd" %in% names(par)){
    if (!is.null(par[["genX.sd"]])) {
      sd.vec[n.peaks] <- par[["genX.sd"]]
    }
  }

  if (!is.null(opt.env)) {
    opt.env[["means"]]   <- rbind(opt.env[["means"]], mean.vec)
    opt.env[["sd"]]      <- rbind(opt.env[["sd"]], sd.vec)
    opt.env[["summits"]] <- rbind(opt.env[["summits"]], summit.vec)
  }

  if (verbose) {
    cat("#------------------------------------------------------\n")
    cat("# New block\n")
    cat("means - wrap: ", mean.vec, "\n")
    cat("sd - wrap: ", sd.vec, "\n")
    cat("summits - wrap: ", summit.vec, "\n")
  }

  if (type == "prolif_model") {
    return(prolif_model(mean.vec, sd.vec, summit.vec, x))

  } else if (type == "density") {
    return(prolif_model_density(mean.vec, sd.vec, summit.vec, x, verbose=verbose, log=log, opt.env=opt.env))

  } else if (type == "neg_log_likelihood") {
    score <- prolif_model_nll(mean.vec, sd.vec, summit.vec, x, verbose=verbose, invert=invert, log=log, opt.env=opt.env)
    if (!is.null(opt.env)) {
      opt.env[["score"]]   <- c(opt.env[["score"]], score)
    }
    return(score)
  } else {

    e <- simpleError("Invalid type argument.")
    stop(e)
  }
}

#-------------------------------------------------------------------------------
#' Residuals between observed CTV trace and Gaussian prolif model
#'
#' @param par named numeric vector of model parameters (see \code{prolif_model_wrapper}).
#' @param n.peaks integer number of peaks in the model.
#' @param x numeric vector of histogram midpoints (log10 scale).
#' @param y numeric vector of observed counts (smoothed or raw) matching x in length.
#' @param fixed named numeric vector of fixed parameters (default NULL).
#' @param opt.env optional environment for logging parameter traces (default NULL).
#' @returns numeric vector of residuals (observed minus predicted) for use with
#'   \code{\link[minpack.lm]{nls.lm}}.
prolif_resid <- function(par, n.peaks, x, y, fixed=NULL, opt.env=NULL) {
  residuals <- y - prolif_model_wrapper(par, n.peaks, x, fixed, type="prolif_model", opt.env=opt.env)
  return(residuals)
}
