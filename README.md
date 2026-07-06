# ProliferationAnalysis

Fits Gaussian mixture models to flow cytometry proliferation traces from division-tracking dyes such as CTV (Cell Trace Violet) and CFSE. Provides automated estimation of peak number and position, followed by optimisation via non-linear least squares (Levenberg-Marquardt, mode="LS", default) or maximum likelihood estimation (mode="MLE"). For the non-linear least sqaures fit, the data is binned and smoothed, for MLE the full likelihood is fit to the raw data. Practically if the fit is good, there is very little difference. The advantage is LS is much faster.  

Supports modelling of a dye-negative cell population (peak X) to account for mixed stained and unstained samples. Calculates standard proliferation statistics including proliferation index, division index, and percent divided.

The input is expected to be raw per-event intensities of the division dye. If you have data on the log scale already, raise it by the exponent first. It is reccomended to pre-gate events to live cells prior to running, but trimming of the CTV+ and CTV- region is generally not needed. For gating CytoExploreR is a good option which allows for templating and processing large volumes of data

# Installation
To install, simply run
```
remotes::install_github("TrynkaLab/ProliferationAnalysis")
```

Code is all in base R except for the Levanberg-Marquad algorithm which uses
the implementation in minpack.lm

# Basic usage

The following fits a proliferation trace using provided data

```
library(ProliferationAnalysis)
data("ctv_example")

# Set the lowest point where you expect undivded cells
gen0.cutoff <- 5.3

plot(density(log10(sample.a)), type="l")
abline(v=gen0.cutoff, lty=2)

# Fit the model
fit <- fit_peaks(trace=sample.a,
                peak.0.lower.bound=gen0.cutoff)
                 

# Get the proliferation statistics
get_prolif_stats(fit)
```


For more details, see the vingette: <a href="https://html-preview.github.io/?url=https://github.com/TrynkaLab/ProliferationAnalysis/blob/main/vignettes/fitting_proliferation_model.html" >Fitting a proliferation model with complex mixture</a>

See `?fit_peaks` and `?get_prolif_stats` for details on all the parameters.