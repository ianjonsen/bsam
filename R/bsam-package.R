

#' Fit Bayesian state-space models to Argos satellite tracking data
#' 
#' Currently models provided are DCRW (for location filtering), DCRWS (for
#' location filtering and behavioural state estimation), and hDCRWS (a
#' hierarchical model for location filtering and behavioural state estimation
#' across multiple animals). The models are fit in JAGS using Markov chain
#' Monte Carlo simulation methods.
#' 
#' \tabular{ll}{ Package: \tab bsam\cr Type: \tab Package\cr Version: \tab
#' 0.43-1\cr Date: \tab 2014-10-07\cr License: \tab GPL-2\cr LazyLoad: \tab
#' yes\cr } Fit Bayesian state-space models to Argos satellite tracking data.
#' Models provided are DCRW - for location filtering; DCRWS - for location
#' filtering and behavioural state estimation with 2 behavioural states; DCRW3S
#' - same as DCRWS but with 3 behavioural states; hDCRW and hDCRWS -
#' hierarchical models for location filtering only, and location filtering with
#' behavioural state estimation, respectively, across multiple animals.
#' 
#' Note that convergence with the DCRW3S is tricky and will generally only
#' converge when there is very clear evidence of more than two movement states
#' in the observed data.
#' 
#' The hierarchical models often provide improved location and/or behavioural
#' state estimates compared to fitting DCRW/DCRWS to individual datasets.
#' 
#' @name bsam-package
#' @aliases bsam-package bsam
#' @docType package
#' @author Ian Jonsen
#' 
#' Maintainer: Ian Jonsen <ian.jonsen@@mq.edu.au>
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
#' @keywords bsam
#' @examples
#' 
#' #data(lbt)
#' #fit = fitSSM(lbt, model="DCRW", tstep=1, adapt=60000, samples=30000, thin=30, chains=2)
#' @importFrom grDevices colorRampPalette dev.off extendrange grey pdf rgb
#' @rawRd 
#' if(.Platform$OS.type == "windows") {
#' importFrom(grDevices, windows)
#' }
#' if(.Platform$OS.type == "unix") {
#' importFrom(grDevices, quartz)
#' }
#' @importFrom graphics axis layout lines matpoints mtext par plot points rug
#' @importFrom stats IQR approx density median quantile rbinom rt runif sd update
#' @importFrom utils data
NULL



##' @name lbt
##' @docType data
##' @title lbt
##' @format csv
##' @keywords data
##' @description lbt data
NULL


