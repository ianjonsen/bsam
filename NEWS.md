# bsam 1.0.1
* removed sp and rworldxtra from Imports

* added BugReports URL to DESCRIPTION

# bsam 1.0.0 

* Simplified movement models by removing the mean turn angle parameter. This tends to improve convergence for the behavioural
switching models.

* Simplified the regularisation / interpolation in the observation models.

* Simplified data preparation code.

* Diagnostic plots (renamed from `diagSSM` to `diag_ssm`) now include the Gelman-Rubin-Brooks shrink factor plots for each parameter.

* New mapping function (`map_ssm`) uses coastline data from `rworldxtra` and `ggplot2` for core plotting functions.

* New plot function (`plot_fit`) to inspect fit to location data.

* Renamed core function `fitSSM` to `fit_ssm`.

* Improved selection of random initial values for MCMC sampling.

* Initial values for location states are now based on a loess smooth through the observed locations. Users can control the degree of smoothing via the `span` argument to `fit_ssm`.


# bsam 0.43.1 (pre-CRAN release)

* ported from source 2016-05-27 mdsumner@gmail.com

* converted to use roxygen2

* Added a `NEWS.md` file to track changes to the package.

```R
f <- "http://web.science.mq.edu.au/~ijonsen/code/bsam_0.43.1.tar.gz"
download.file(f, basename(f), mode = "wb")
system(sprintf("tar zxvf %s", basename(f)))
Rd2roxygen::Rd2roxygen("bsam")
```







