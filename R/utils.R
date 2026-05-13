#-------------------------------------------------------------------------------
#' Find Local Max
#'
#' Find the maximum value of y within a range of x values
#'
#' @param x vector of numbers matching y in length.
#' @param y vector of numbers matching x in length.
#' @param lower.x lower limit of x to search.
#' @param upper.x upper limit of x to search.
#' @returns the value of x where y is maximised.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' find_local_max(x, y, -1, 1)
find_local_max <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymax <- max(y.f)
  return(mean(x.f[y.f==ymax]))
}

#-------------------------------------------------------------------------------
#' Find Local Min
#'
#' Find the minimum value y of x within a range of x values
#'
#' @param x vector of numbers matching y in length.
#' @param y vector of numbers matching x in length.
#' @param lower.x lower limit of x to search.
#' @param upper.x upper limit of x to search.
#' @returns the value of x where y is minimized
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' find_local_min(x, y, -1, 1)
find_local_min <- function(x, y, lower.x, upper.x){
  y.f <- y[x >= lower.x & x <= upper.x]
  x.f <- x[x >= lower.x & x <= upper.x]
  ymin <- min(y.f)
  return(mean(x.f[y.f==ymin]))
}

#-------------------------------------------------------------------------------
#' Nearest Neighbour smoother
#'
#' Smooth a trace using a kernel of size window.
#'
#' @param y a vector y
#' @param window how many values up and down are averaged
#' @returns the mean of each y value and the surrounding values specified by window
#'
#' @examples
#' y <- rnorm(100)
#' nn_smoother(y, 2)
#' @export
nn_smoother <- function(y, window=2) {
  y.out <- c()

  for(i in 1:length(y)) {
    i.min <- i-window
    i.min <- ifelse(i.min < 0, 0, i.min)
    i.max <- i+window
    i.max <- ifelse(i.max > length(y), length(y), i.max)

    y.out[i] <- mean(y[i.min:i.max])
  }
  return(y.out)
}

#-------------------------------------------------------------------------------
#' Find nearest index
#'
#' @description
#' Find the index of the value in x that is closest to value
#'
#' @param x a vector
#' @param value the value to search the nearest value in x
#' @returns the index in x closest to value
#'
#' @examples
#' x <- rnorm(100)
#' nearest_index(x, 0)
nearest_index <- function(x, value) {
  which.min(abs(x-value))
}
