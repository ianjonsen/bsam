# bsam 1.1.2
* fixed indexing on priors for initial location state. A sufficiently short time step could cause the hierarchical models to crash on compilation.

# bsam 1.1.1
* fixed the mis-ordering of animal id's in the summary and data output data.frames caused when a hierarchical model was fit to a dataset with >= 10 individuals

# bsam 1.1.0.9000
* changed diag_ssm plot so that diagnostic panels for psi parameters are split into multiple pages of
5 rows/page when fitting a hierarchical model. Previous version packed all the psi diagnostic plots into a
single page

# bsam 1.1.0
* major bug fix - data indexing when fitting hierarchical models was causing improper fits to multi-individual datasets with > 2 individuals. Indexing now works as intended 
* Added simulate function 
* Prevented individual animal tracks from being re-ordered by split function so output track order now matches input track order

# bsam 1.0.1.9003
* added list of initial values to the output list object. Making these available can aid diagnosing lack of convergence and inform choice of span argument value for generating location state initial values

# bsam 1.0.1.9002
* added `get_summary` function to extract summary `data_frame` from `fit_ssm` output objects. The `data_frame` can optionally be written to a .csv file

# bsam 1.0.1.9001
* added `sp` and `rworldxtra` back to Imports list in DESCRIPTION. Ensures all required packages are installed, otherwise `map_ssm` will return an error when attempting to load `countriesHigh` data if `sp` and/or `rworldxtra` are not installed

# bsam 1.0.1.9000
* removed sp and rworldxtra from Imports

* added BugReports URL to DESCRIPTION

# bsam 1.0.0 

* Simplified movement models by removing the mean turn angle parameter. This tends to improve convergence for the behavioural switching models

* Simplified the regularisation / interpolation in the observation models

* Simplified data preparation code

* Diagnostic plots (renamed from `diagSSM` to `diag_ssm`) now include the Gelman-Rubin-Brooks shrink factor plots for each parameter

* New mapping function (`map_ssm`) uses coastline data from `rworldxtra` and `ggplot2` for core plotting functions

* New plot function (`plot_fit`) to inspect fit to location data

* Renamed core function `fitSSM` to `fit_ssm`

* Improved selection of random initial values for MCMC sampling

* Initial values for location states are now based on a loess smooth through the observed locations. Users can control the degree of smoothing via the `span` argument to `fit_ssm`


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







