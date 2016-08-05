# bsam

**bsam** - **B**ayesian **S**tate-space models for **A**nimal **M**ovement. 

bsam is an R package that fits Bayesian state-space models via JAGS to Argos satellite tracking data. The models filter error-prone Argos locations and estimate behavioural states associated with two fundamentally different movement patterns (directed, fast movements and relatively undirected, slow movements). The models can be fit to individual animal tracks or simultaneously to a group of individuals. Plotting functions are provided to help assess lack of MCMC convergence, map estimated tracks and visualise fit to the observations.

Read `?bsam` for more details on functionality and `?fit_ssm` for details and examples of how to use the package. 

## Installation

First ensure that you have a working copy of JAGS (>= 4.2.0) for the [rjags](https://cran.r-project.org/package=rjags) package, 
see instructions below. 

```R
library(rjags)
```

Get the released version of bsam from CRAN:

```R
install.packages("bsam")
```

Or the development version from GitHub:

```R
# install.packages("devtools")  
devtools::install_github("ianjonsen/bsam")
```

### JAGS installation

Install JAGS: http://mcmc-jags.sourceforge.net/

