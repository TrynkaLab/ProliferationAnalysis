# Deprecated aliases for renamed public API functions.
# These wrappers emit a deprecation warning and forward to the new names.

#' @rdname nn_smoother
#' @export
nn.smoother <- function(...) {
  .Deprecated("nn_smoother")
  nn_smoother(...)
}

#' @rdname fit_peaks
#' @export
fit.peaks <- function(...) {
  .Deprecated("fit_peaks")
  fit_peaks(...)
}

#' @rdname find_peak_x_approx_mode
#' @export
find.peak.x.approx.mode <- function(...) {
  .Deprecated("find_peak_x_approx_mode")
  find_peak_x_approx_mode(...)
}

#' @rdname get_prolif_stats
#' @export
get.prolif.stats <- function(...) {
  .Deprecated("get_prolif_stats")
  get_prolif_stats(...)
}
