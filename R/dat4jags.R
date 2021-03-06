##' Format track data for filtering
##'
##' This is an internal function used by \code{fit_ssm} to format track
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
##' @param tpar generalised t-distribution parameters for ARGOS location classes. By 
##' default dat4jags uses the parameters estimated in Jonsen et al (2005) Ecology 86:2874-2880
##' but users may specify other ARGOS error parameter values via the \code{tpar} function.
##' @return A list with components
##' \item{\code{id}}{the unique identifier for each dataset}
##' \item{\code{y}}{a 2 column matrix of the lon,lat observations}
##' \item{\code{itau2}}{a 2 column matrix of the ARGOS precision (1/scale) parameters}
##' \item{\code{nu}}{a 2 column matrix of the ARGOS df parameters}
##' \item{\code{idx}}{a vector of interpolation indices}
##' \item{\code{ws}}{a vector of interpolation weights}
##' \item{\code{ts}}{the times at which states are predicted (POSIXct,GMT)}
##' \item{\code{obs}}{the input observed data frame}
##' \item{\code{tstep}}{the time step specified in the \code{fitSSM} call}
##' @references Jonsen ID, Mills Flemming J, Myers RA (2005) Robust state-space modeling of
##' animal movement data. Ecology 86:2874-2880 (Appendix A)
##' @importFrom dplyr left_join
##' @importFrom lubridate ymd_hms
##' @export
dat4jags <- function (d, tstep=1, tpar=tpar()) {
  
  ## Check input data columns; return error if unexpected number
  if(ncol(d) != 5 && ncol(d) != 7) {
    stop("Input data should have 5 or 7 columns, fitSSM help for details")
  }
  if(ncol(d) == 7 && (names(d)[6] != "lonerr" || names(d)[7] != "laterr")) {
    stop("Input data columns 6 and 7 must be labelled `lonerr` and `laterr`, respectively")
  }
  
  ## Ensure POSIXct dates
  d$date <- ymd_hms(d$date, tz = "GMT")

  ## Check ARGOS lc - make sure it has class 'character' prior to left_join w tpar
  d$lc <- as.character(d$lc)
  
  ## Merge ARGOS error (t-distribution) fixed parameters
  dnew <- left_join(d, tpar, by = "lc")
  
  if(ncol(d) == 7) {
    dnew$itau2.lon <- d$lonerr ^ -2
    dnew$itau2.lat <- d$laterr ^ -2
    dnew <- dnew[, -c(6,7)]
  }
  
  ## now coerce lc to an ordered factor
  dnew$lc <- factor(as.character(dnew$lc), levels=c("3", "2", "1", "0", "A", "B", 
                                              "Z", "F", "G"), 
                 ordered=TRUE)
  
  dostuff <- function(x) {
    ## Interpolation indices and weights
    dt <- tstep * 86400
    tms <- (as.numeric(x$date) - as.numeric(x$date[1])) / dt
    index <- floor(tms) + 1
    weights <- 1 - (tms - (index - 1))
  
    list(
      id = x$id[1],
      y = cbind(x$lon, x$lat),
      itau2 = cbind(x$itau2.lon, x$itau2.lat),
      nu = cbind(x$nu.lon, x$nu.lat),
      idx = index,
      ws = weights,
      ts = seq(x$date[1], by = dt, length.out = max(index)),
      obs = x,
      tstep = tstep
    )
  }
  
  lapply(split(dnew, dnew$id), dostuff)
}