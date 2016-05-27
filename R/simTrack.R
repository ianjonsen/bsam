#' @param theta 
#'
#' @param gamma 
#' @param alpha 
#' @param vcov 
#'
#' @export
`simTrack` = function(T = 100, theta = c(0, pi), gamma = c(0.9, 0.2), alpha = c(0.2, 0.8), vcov = c(0.01, 0, 0, 0.02)){

start.date <- strptime(format(Sys.time(), "%d/%m/%y %H:%M:%S"), "%d/%m/%y %H:%M:%S", tz = "GMT")
x <- matrix(NA, T, 2)
Tdx <- matrix(NA, T-1, 2)
x.mn <- matrix(NA, T-1, 2)
b <- c()
mu <- c()
tau.lon <- c()
tau.lat <- c()
nu.lon <- c()
nu.lat <- c()
Sigma <- matrix(vcov, 2, 2)

## Start as random walk
x[1, ] <- mvrnorm(1, c(-65,40), Sigma) #randomize starting position
x[2, ] <- mvrnorm(1, x[1,], Sigma)

## set first behav state
b[1] <- rbinom(1, 1, 0.5) + 1

## NOTES:
##
##  It may be easier to start with simulating the switch from state 1 to state 2 as a fn of environment. This may be a better reflection of reality, as turtles should use other (potentially internal) cues for leaving the slow movement state.

## iterate DCRWS for successive steps
for(i in 2:(T-1)){
	b[i] = rbinom(1, 1, alpha[b[i-1]]) + 1
	
	Tdx[i,1] <- cos(theta[b[i]]) * (x[i,1] - x[i-1,1]) + sin(theta[b[i]]) * (x[i,2] - x[i-1,2])
	Tdx[i,2] <- -sin(theta[b[i]]) * (x[i,1] - x[i-1,1]) + cos(theta[b[i]]) * (x[i,2] - x[i-1,2])
	
	x.mn[i,] <- x[i,] + Tdx[i,] * gamma[b[i]]
	x[i+1,] <- mvrnorm(1, x.mn[i,], Sigma)
	}

## Behav state at last time step
b[T] = rbinom(1, 1, alpha[b[T-1]]) + 1

## Randomly draw lc's, using class proportions from LBT data as the probability vector
lc = factor(sample(c(3,2,1,0,"A","B"), T, replace=TRUE, prob=c(0.03, 0.09, 0.16, 0.14, 0.26, 0.32)), levels=c(3,2,1,0,"A","B"), ordered=TRUE)
	
## Add measurement error to true locations
## Error SD's and df's from Jonsen et al. (2005) Ecology
tau.lon[1] <- 0.2898660 / 6366.71 * 180 / pi
tau.lon[2] <- 0.3119293 / 6366.71 * 180 / pi
tau.lon[3] <- 0.9020423 / 6366.71 * 180 / pi
tau.lon[4] <- 2.1625936 / 6366.71 * 180 / pi
tau.lon[5] <- 0.5072920 / 6366.71 * 180 / pi
tau.lon[6] <- 4.2050261 / 6366.71 * 180 / pi

tau.lat[1] <- 0.1220553 / 6366.71 * 180 / pi
tau.lat[2] <- 0.2605126 / 6366.71 * 180 / pi
tau.lat[3] <- 0.4603374 / 6366.71 * 180 / pi
tau.lat[4] <- 1.607056 / 6366.71 * 180 / pi
tau.lat[5] <- 0.5105468 / 6366.71 * 180 / pi
tau.lat[6] <- 3.041276 / 6366.71 * 180 / pi

nu.lon[1] <- 3.070609
nu.lon[2] <- 1.220822
nu.lon[3] <- 2.298819
nu.lon[4] <- 0.9136517
nu.lon[5] <- 0.786954
nu.lon[6] <- 1.079216

nu.lat[1] <- 2.075642
nu.lat[2] <- 6.314726
nu.lat[3] <- 3.896554
nu.lat[4] <- 1.010729
nu.lat[5] <- 1.057779
nu.lat[6] <- 1.331283

lon.obs = x[,1] + tau.lon[as.numeric(lc)] * rt(T, nu.lon[as.numeric(lc)])
lat.obs = x[,2] + tau.lat[as.numeric(lc)] * rt(T, nu.lat[as.numeric(lc)])

## time interval is 10 minutes, or 600 s
dates = seq(start.date, start.date + (T-1) * 600, by=600)

simdat = data.frame(date=dates, lon=x[,1], lat = x[,2], lon.obs, lat.obs, lc, b)

subdat = simdat[sort(sample(1:T, T/4, replace=FALSE)), ]

ssmdat = with(subdat, data.frame(id=rep(1, nrow(subdat)), date=as.character(date), lc, lon=lon.obs, lat=lat.obs))


# return simulated data and SSM results
list(simdat=simdat, ssmdat=ssmdat)
}