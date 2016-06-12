#' Fit Bayesian state-space models to animal tracking data
#' 
#' Models provided are DCRW (for location filtering), DCRWS (for
#' location filtering and behavioural state estimation), and their hierarchical 
#' versions (hDCRW, hDCRWS) to estimate parameters jointly across multiple 
#' individual tracking datasets. The models are fit in JAGS using Markov chain
#' Monte Carlo simulation methods. The models are intended to be fit to Argos
#' satellite tracking data but options exist to allow fits to other tracking
#' data types (type \code{?fit_ssm} for details).
#' 
#' \tabular{ll}{ Package: \tab bsam\cr Type: \tab Package\cr Version: \tab
#' 0.50.0\cr Date: \tab 2016-06-09\cr License: \tab GPL-2\cr LazyLoad: \tab
#' yes\cr } Fit Bayesian state-space models to Argos satellite tracking data.
#' Models provided are DCRW - for location filtering; DCRWS - for location
#' filtering and behavioural state estimation with 2 behavioural states; hDCRW 
#' and hDCRWS - hierarchical models for location filtering only, and location 
#' filtering with behavioural state estimation, respectively, across multiple 
#' animals.
#' 
#' The hierarchical models can provide improved location and/or behavioural
#' state estimates compared to fitting DCRW/DCRWS to individual datasets.
#' 
#' @name bsam-package
#' @aliases bsam-package bsam
#' @docType package
#' @author Ian Jonsen
#' 
#' Maintainer: Ian Jonsen <ian.jonsen@mq.edu.au>
#' @seealso fitSSM
#' @references Jonsen ID, Myers RA, Mills Flemming J (2003) Meta-analysis of
#' animal movement using state-space models. Ecology 84:3055-3063
#' 
#' Jonsen ID, Mills Flemming J, Myers RA (2005) Robust state-space modeling of
#' animal movement data. Ecology 86:2874-2880
#' 
#' Jonsen ID, Myers RA, James MC (2007) Identifying leatherback turtle foraging
#' behaviour from satellite telemetry using a switching state-space model.
#' Marine Ecology Progress Series 337:255-263
#' 
#' Jonsen ID (2016) Joint estimation over multiple individuals improves behavioural 
#' state inference from animal movement data. Scientific Reports 6:20625
#' @keywords bsam
#' @examples
#' \dontrun{
#' # Fit DCRW model for state filtering and regularization
#' data(ellie)
#' fit <- fit_ssm(ellie, model = "DCRW", tstep = 1, adapt = 5000, samples = 5000, 
#'               thin = 5, span=0.2)
#' diag_ssm(fit)
#' map_ssm(fit)
#' plot_fit(fit)
#' }
#' @importFrom graphics axis layout lines matpoints mtext par plot points rug
#' @importFrom stats IQR approx density median quantile rbinom rt runif sd update
#' @importFrom utils data
#' @importFrom grDevices extendrange grey
#' @importFrom stats cov loess loess.control predict rbeta rlnorm rnorm
#' @importFrom gridExtra grid.arrange
NULL



##' @name ellie
##' @docType data
##' @title Elephant seal Argos satellite data (2 individuals)
##' @format .RData
##' @keywords data 
##' @description Example elephant seal Argos tracking data. Data were sourced from 
##' the Integrated Marine Observing System (IMOS) - IMOS is supported by the 
##' Australian Government through the National Collaborative Research Infrastructure 
##' Strategy and the Super Science Initiative.
NULL


