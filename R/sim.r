#' Simulate from the DCRW model with Argos location errors
#'
#' For testing bsam models
#'
#' @param Nt number of time steps to simulate
#' @param gamma move persistence parameter
#' @param Sigma variance-covariance matrix for movement process
#' @return a data_frame of true locations and locations with Argos error
#' @importFrom mvtnorm rmvnorm
#' @importFrom tibble data_frame
#' @export
sim <- function(Nt = 100,
                gamma = 0.8,
                Sigma = matrix(c(0.01, 0, 0, 0.01), 2, 2),
                amf = tpar()) {
  
  start.date = Sys.time()
  
  Y <- X <- matrix(NA, Nt, 2)
  
  ## Start as random walk
  X[1,] <- c(0, 0) 
  X[2,] <- rmvnorm(1, X[1, ], Sigma)
  
  ## iterate DCRWS for successive steps
  for (i in 2:(Nt - 1)) {
    X[i + 1, ] = rmvnorm(1, X[i, ] + (X[i, ] - X[i - 1, ]) * gamma, Sigma)
  }
  
  ## time interval is nominally 1 h
  dt <- seq(start.date, start.date + (Nt - 1) * 1/24 * 86400,
            by = 1/24 * 86400)
   
  ## probability vector
  lc <- factor(
    sample(
      c(3, 2, 1, 0, "A", "B"),
      Nt,
      replace = TRUE,
      prob = c(0.03, 0.04, 0.059, 0.145, 0.371, 0.353)
    ),
    levels = c(3, 2, 1, 0, "A", "B"),
    ordered = TRUE
  )

  d <- data_frame(
    date = dt,
    lon = X[, 1],
    lat = X[, 2],
    lc = lc
  )
  
  ## Merge ARGOS error multiplication factors
  d <- merge(d, amf, by = "lc", all.x = TRUE)
  d <- d[order(d$date), ]
  
  lo <- d$lon + sapply(1:Nt, function(i) {
    rt(1, d$nu.lon[i]) * d$itau2.lon[i]^-0.5
  })
  la <- d$lat + sapply(1:Nt, function(i) {
    rt(1, d$nu.lat[i]) * d$itau2.lat[i]^-0.5
  })
  
  with(d, data_frame(id = round(runif(1, 1, 10000)), date, lc, lon.o = lo, lat.o = la, lon, lat))
  
}
