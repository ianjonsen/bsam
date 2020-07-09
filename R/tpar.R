##' ARGOS Error Fixed Parameters for Location Classes
##'
##' This is an internal function used by \code{dat4jags} to specify measurement
##'   error parameters.
##' 
##' These are the fixed parameters (t-distribution scale & df) for ARGOS
##' error classes, from Jonsen et al (2005) Ecology 86:2874-2880. 
##' @title ARGOS Error Fixed Parameters
##' @return A dataframe with columns
##' \item{\code{lc}}{ARGOS location class as an ordered factor}
##' \item{\code{itau2.lon}}{precision parameters for longitude in degrees}
##' \item{\code{itau2.lat}}{precision parameters for latitude in degrees}
##' \item{\code{nu.lon}}{df parameters for longitude}
##' \item{\code{nu.lat}}{df parameters for latitude}
##' @export
tpar <- function() {
  data.frame(lc = c("3", "2", "1", "0", "A", "B", "Z", "F", "G"),
             itau2.lon = (c(0.289866, 0.3119293, 0.9020423, 2.1625936, 0.507292, 
                       4.2050261, 4.2050261, 0.01, NA) / 6366.71 * 180 / pi) ^ -2,
             itau2.lat = (c(0.1220553, 0.2605126, 0.4603374, 1.607056, 0.5105468, 
                       3.041276, 3.041276, 0.01, NA) / 6366.71 * 180 / pi) ^ -2,
             nu.lon = c(3.070609, 1.220822, 2.298819, 0.9136517, 0.786954, 
                      1.079216, 1.079216, 100000, 100000),
             nu.lat = c(2.075642, 6.314726, 3.896554, 1.010729, 1.057779, 
                      1.331283, 1.331283, 100000, 100000))
}
