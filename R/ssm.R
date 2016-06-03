#' Fits state-space model to Argos data
#' 
#' Takes output from dat4jags, sets up initial values, calls JAGS, and
#' aggregates results. Intended for internal use, called by \code{fitSSM}.
#' 
#' @param d structured data from dat4jags to be passed to JAGS
#' @param model the state-space model to be fit: DCRW or DCRWS
#' @param adapt number of samples in adaptation/burnin phase 
#' @param samples number of posterior samples
#' @param thin thinning factor to reduce posterior sample autocorrelation
#' @param chains number of parallel McMC chains to run
#' the only model options.
#' @return Returns a list of McMC samples from marginal posteriors and a
#' summary data.frame of mean and median position estimates.
#' @seealso Function to be called by \code{\link{fitSSM}}.
#' @importFrom rjags jags.samples
#' @export
ssm <- function (d, model = "DCRW", adapt, samples, thin, chains, span)
{
    ssm1 = function(dd) {
      gamma <- 0.5
      fit.lon <- loess(lon ~ as.numeric(date), data=dd$obs, span=span,
                       na.action="na.exclude", control=loess.control(surface="direct"))
      fit.lat <- loess(lat ~ as.numeric(date), data=dd$obs, span=span,
                       na.action="na.exclude", control=loess.control(surface="direct"))
      ## Predict track, increments and stochastic innovations
      xs <- cbind(predict(fit.lon, newdata=data.frame(date=as.numeric(dd$ts))),
                  predict(fit.lat, newdata=data.frame(date=as.numeric(dd$ts))))
      ds <- xs[-1, ] - xs[-nrow(xs), ]
      es <- ds[-1, ] - gamma*ds[-nrow(ds), ]
      
      ## Estimate process variance components for initial values
      V <- cov(es)
      isigma2 <- diag(V)^-1
      rho <- V[1,2] / prod(sqrt(isigma2))
      iSigma <- matrix(c(isigma2[1], rho, rho, isigma2[2]), 2, 2)
      
      data <- with(dd, list(y = y, idx = idx, w = ws, itau2 = itau2, nu = nu, 
                            Nx = nrow(xs), Ny=nrow(dd$y)))
      ## inits
      switch(model,
        	DCRW = {
            inits <- list(list(iSigma = iSigma, gamma = 0.5,
               theta = 0, x = xs, logpsi = runif(1,-1,1)),
              list(iSigma = iSigma, gamma = 0.4, theta = 0.1, x = xs,
               	logpsi = runif(1,-1,1)), list(iSigma = iSigma, gamma = 0.45,
               	theta = -0.1, x = xs, logpsi = runif(1,-1,1)),
              list(iSigma = iSigma, gamma = 0.3, theta = -0.2, x = xs,
               	logpsi = runif(1,-1,1))) 
            params <- c("Sigma", "x", "theta", "gamma", "psi")
        		},
        	DCRWS = {
        	  b <- rep(1, nrow(xs))
        	  d <- sqrt(ds[,1]^2+ds[,2]^2)
        	  b[d > median(d)] <- 2
            inits <- list(list(iSigma = iSigma, gamma = c(0.8, NA), dev = 0.6,
            	  tmp = c(0.45, 0.55), alpha = c(0.55, 0.45), lambda = c(0.45, NA),
            	  b = b, x = xs, logpsi = runif(1,-1,1)),
              list(iSigma = iSigma, gamma = c(0.85, NA), dev = 0.65,
            	  tmp = c(0.6, 0.4), alpha = c(0.4, 0.6), lambda = c(0.5, NA),
            	  b = b, x = xs, logpsi = runif(1,-1,1)), 
            	list(iSigma = iSigma, gamma = c(0.75, NA), dev = 0.7,
            	  tmp = c(0.5, 0.6), alpha = c(0.5, 0.5), lambda = c(0.55, NA),
            	  b = b, x = xs, logpsi = runif(1,-1,1)),
            	list(iSigma = iSigma, gamma = c(0.6, NA), dev = 0.5,
            	  tmp = c(0.6, 0.5), alpha = c(0.65, 0.45), lambda = c(0.4, NA),
            	  b = b, x = xs, logpsi = runif(1,-1,1)))
            params <- c("Sigma", "x", "theta", "gamma", "alpha", "b", "psi")
        		})

	if(chains==1) inits <- inits[[1]]
	if(chains==2) inits <- list(inits[[1]],inits[[2]])
	if(chains==3) inits <- list(inits[[1]],inits[[2]],inits[[3]])

	model.file <- file.path(system.file("jags", package="bsam"), paste(model, ".txt", sep=""))

	burn <- rjags::jags.model(model.file, data, inits, n.chains=chains, n.adapt=adapt/2)
	update(burn, n.iter=adapt/2)
	psamples <- rjags::jags.samples(burn, params, n.iter=samples, thin=thin)

	lon <- apply(psamples$x[,1,,],1, mean)
	lat <- apply(psamples$x[,2,,],1, mean)
	lon.q <- apply(psamples$x[,1,,],1, quantile, c(0.025, 0.5, 0.975))
	lat.q <- apply(psamples$x[,2,,],1, quantile, c(0.025, 0.5, 0.975))

	summary <- data.frame(id = dd$id, date = dd$ts,
	                     lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
	                     lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,])
	model <- model
	mcmc.settings <- list(burnin = adapt, posterior.samples = samples, thinning = thin, n.chains = chains)
	data <- dd$obs

	if(model == "DCRWS") {
	  b = apply(psamples$b, 1, mean)
	  b.5 = apply(psamples$b, 1, median)
	  summary <- data.frame(summary, b=b, b.5=b.5)
	}

	out <- list(summary=summary, mcmc=psamples, model=model, mcmc.settings=mcmc.settings,
		timestep=dd$tstep, Nx=nrow(xs), data=data)
	out
    }
    lapply(d, ssm1)
}
