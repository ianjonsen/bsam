##' Format track data for filtering
##'
##' This is an internal function used by \code{fitSSM} to format track
##' data for JAGS.
##'
##' The input track is given as a dataframe where each row is an
##' observed location and columns
##' \describe{
##' \item{'id'}{individual animal identifier,}
##' \item{'date'}{observation time (POSIXct,GMT),}
##' \item{'lc'}{ARGOS location class,}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude.}
##' }
##' Location classes can include Z, F, and G; where the latter two 
##' are used to designate fixed (known) locations (e.g. GPS locations)
##' and "generic" locations (e.g. geolocation data) where the user 
##' supplies the error standard deviations, either via the 
##' \code{tpar} function or as two extra columns in the input data.
##' 
##' From this \code{dat4jags} calculates interpolation indices \code{idx} and
##' weights \code{ws} such that if \code{x} is the matrix of predicted
##' states, the fitted locations are \code{ws*x[idx+1,] +
##' (1-ws)*x[idx+2,]}. 
##'
##' @title Correlated Random Walk Filter
##' @param d a data frame of observations (see details)
##' @param tstep the time step to predict to (in days)
##' @param extrap if TRUE, the final predicted state occurs
##'   immediately before the last observation, otherwise the final
##'   predicted state occurs immediately after the last observation.
##' @param tpar generalised t-distribution parameters for ARGOS location classes
##' @return A list with components
##' \item{\code{y}}{a 2 column matrix of the lon,lat observations}
##' \item{\code{K}}{a 2 column matrix of the ARGOS scale factors}
##' \item{\code{idx}}{a vector of interpolation indices}
##' \item{\code{ws}}{a vector of interpolation weights}
##' \item{\code{ts}}{the times at which states are predicted (POSIXct,GMT)}
##' \item{\code{dt}}{the time step at which states are predicted (secs)}
##' @references Jonsen ID, Mills Flemming J, Myers RA (2005) Robust state-space modeling of
#' animal movement data. Ecology 86:2874-2880 (Appendix A)
##' @export
dat4jags <- function (d, tstep=1, tpar=tpar()) {
  
  ## Check ARGOS location accuracies
  d$lc <- factor(as.character(d$lc), levels=c("3", "2", "1", "0", "A", "B", 
                                              "Z"), 
                 ordered=TRUE)
  ## Ensure POSIXct dates
  d$date <- as.POSIXct(d$date, format="%Y-%m-%d %H:%M:%S", tz="GMT")
    
  ## Merge ARGOS error (t-distribution) fixed parameters
  d <- merge(d, tpar, by="lc", all.x=TRUE)
  d <- d[order(d$date), ]
  
  dostuff <- function(dd) {
    ## Interpolation indices and weights
    dt <- tstep * 86400
    tms <- (as.numeric(dd$date) - as.numeric(dd$date[1])) / dt
    index <- floor(tms) + 1
    weights <- 1 - (tms - (index - 1))

    list(id = dd$id[1], y = cbind(dd$lon, dd$lat),
       itau2 = cbind(dd$itau2.lon, dd$itau2.lat),
       nu = cbind(dd$nu.lon, dd$nu.lat),
       idx = index,
       ws = weights,
       ts = seq(dd$date[1], by = dt, length.out = max(index)),
       dt = dt, obs = dd, tstep = tstep)
  }
  by(d, d$id, dostuff)
}