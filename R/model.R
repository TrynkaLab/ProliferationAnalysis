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
#' Takes a vector of parameters par, number of peaks, and data X
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
prolif_resid <- function(par, n.peaks, x, y, fixed=NULL, opt.env=NULL) {
  residuals <- y - prolif_model_wrapper(par, n.peaks, x, fixed, type="prolif_model", opt.env=opt.env)
  return(residuals)
}
