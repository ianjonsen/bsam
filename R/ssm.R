#' Fits state-space model to Argos data
#' 
#' Takes output from dat4jags, sets up initial values, calls JAGS, and
#' aggregates results. Intended for internal use, called by fitSSM.
#' 
#' 
#' @param loc.list Data from dat4jags to be passed to JAGS
#' @param model The state-space model to be fit. Currently, DCRW and DCRWS are
#' the only model options.
#' @param \dots Other arguments can be passed.
#' @return Returns a list of McMC samples from marginal posteriors and a
#' summary data.frame of mean and median position estimates.
#' @seealso Function to be called by \code{\link{fitSSM}}.
#' @importFrom rjags jags.samples
#' @export
`ssm` =
function (loc.list, model = "DCRW", adapt, samples, thin, chains, ...)
{
    if (!model %in% c("DCRW","DCRWS","DCRW3S"))
        stop("model not implemented")
    ssm1 = function(input) {
        y = input$y
        idx = input$idx
        RegN = input$RegN
        j = input$j
        x = rbind(y[idx[-length(idx)], ], y[max(idx) - 1, ])
        id = input$id
        row.na = which(is.na(x[, 1]))
        x[row.na, 1] = approx(seq(nrow(x)), x[, 1], xout = row.na, rule = 2)$y
        row.na = which(is.na(x[, 2]))
        x[row.na, 2] = approx(seq(nrow(x)), x[, 2], xout = row.na, rule = 2)$y
        row.na = which(is.na(y[, 1]) | is.na(y[, 2]))
        if (length(row.na) == 0) {
            y.init = y
            y.init[!is.na(y.init)] = NA
        	}
        else {
            y.init = y
            y.init[!is.na(y.init)] = 9999
            yinit.xna = approx(seq(nrow(y)), y[, 1], xout = row.na, rule = 2)$y
            yinit.yna = approx(seq(nrow(y)), y[, 2], xout = row.na, rule = 2)$y
            y.init[row.na, ] = cbind(yinit.xna, yinit.yna)
            y.init[y.init == 9999] = NA
        	}
        start.date = input$first.date
        tstep.sec = input$tstep * 86400
        steplims = seq(start.date, by = paste(tstep.sec, "sec"), length = length(idx))
     	jags.data = list(y = y, itau2 = input$itau2, nu = input$nu, j = j, idx = idx,
     		RegN = RegN)
        iSigma = matrix(c(1, 0, 0, 1), 2, 2)
        switch(model,
        	DCRW = {
            	jags.inits <- list(list(iSigma = iSigma, gamma = 0.5,
               	 theta = 0, x = x, logpsi = runif(1,-1,1), y = y.init),
               	 list(iSigma = iSigma, gamma = 0.4, theta = 0.1, x = x,
               	 logpsi = runif(1,-1,1), y = y.init), list(iSigma = iSigma, gamma = 0.45,
               	 theta = -0.1, x = x, logpsi = runif(1,-1,1), y = y.init),
               	 list(iSigma = iSigma, gamma = 0.3, theta = -0.2, x = x,
               	 logpsi = runif(1,-1,1), y = y.init))
            	jags.params <- c("Sigma", "x", "theta", "gamma", "psi")
        		},
        	DCRWS = {
            	jags.inits <- list(list(iSigma = iSigma, gamma = c(0.8, NA), dev = 0.6,
            	tmp = c(0.45, 0.55), alpha = c(0.55, 0.45), lambda = c(0.45, NA),
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.85, NA), dev = 0.65,
            	tmp = c(0.6, 0.4), alpha = c(0.4, 0.6), lambda = c(0.5, NA),
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.75, NA), dev = 0.7,
            	tmp = c(0.5, 0.6), alpha = c(0.5, 0.5), lambda = c(0.55, NA),
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.6, NA), dev = 0.5,
            	tmp = c(0.6, 0.5), alpha = c(0.65, 0.45), lambda = c(0.4, NA),
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init))
            	jags.params <- c("Sigma", "x", "theta", "gamma", "alpha", "b", "psi")
        		},
        	DCRW3S = {
            	jags.inits <- list(list(iSigma = iSigma, gamma = c(0.8, NA, NA),
            	dev = c(0.6, 0.2), tmp = c(0.45, 0.55, 0.5), alpha1 = c(0.55, 0.45, 0.35),
            	alpha2 = c(0.6, 0.4, 0.2), lambda = c(0.45, NA, NA), lambda1=0.5,
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.85, NA, NA),
            	dev = c(0.65,0.45), tmp = c(0.6, 0.4, 0.35), alpha1 = c(0.4, 0.6, 0.8),
            	alpha2 = c(0.5, 0.5, 0.5), lambda = c(0.5, NA, NA), lambda1 = 0.6,
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.75, NA, NA),
            	dev = c(0.7,0.3), tmp = c(0.5, 0.6, 0.7), alpha1 = c(0.5, 0.5, 0.5),
            	alpha2 = c(0.55, 0.45, 0.35), lambda = c(0.55, NA, NA), lambda1 = 0.4,
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init), list(iSigma = iSigma, gamma = c(0.6, NA, NA),
            	dev = c(0.3, 0.4), tmp = c(0.6, 0.5, 0.3), alpha1 = c(0.65, 0.45, 0.25),
            	alpha2 = c(0.65, 0.45, 0.25), lambda = c(0.4, NA, NA), lambda1 = 0.55,
            	b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(1,-1,1),
            	y = y.init))
            	jags.params <- c("Sigma", "x", "theta", "gamma", "alpha1", "alpha2", "b",
            	"psi")
            })

	if(chains==1) jags.inits = jags.inits[[1]]
	if(chains==2) jags.inits = list(jags.inits[[1]],jags.inits[[2]])
	if(chains==3) jags.inits = list(jags.inits[[1]],jags.inits[[2]],jags.inits[[3]])

	model.file = paste(system.file('jags',package='bsam'), "/", model, ".txt", sep="")

	burn = jags.model(model.file, jags.data, jags.inits, n.chains=chains, n.adapt=adapt/2)
	update(burn, n.iter=adapt/2)
	psamples = rjags::jags.samples(burn, jags.params, n.iter=samples, thin=thin)

	lon = apply(psamples$x[,1,,],1, mean)
	lat = apply(psamples$x[,2,,],1, mean)
	lon.q = apply(psamples$x[,1,,],1, quantile, c(0.025, 0.5, 0.975))
	lat.q = apply(psamples$x[,2,,],1, quantile, c(0.025, 0.5, 0.975))

	switch(model,
		DCRW = {
			summary = data.frame(id=as.character(id), date =
				as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
				lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
				lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,])
			model = model
			mcmc.settings = list(burnin = adapt, posterior.samples = samples, thinning = thin, n.chains = chains)
			step = with(summary, difftime(date[2], date[1],
				units="hours"))
			y = data.frame(lon=y[,1], lat=y[,2])
			},
		DCRWS = {
			b = apply(psamples$b, 1, mean)
			b.5 = apply(psamples$b, 1, median)
			summary = data.frame(id=as.character(id), date =
				as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
				lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
				lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,], b, b.5)
			model = model
			mcmc.settings = list(burnin = adapt, posterior.samples = samples, thinning = thin, n.chains = chains)
			step = with(summary, difftime(date[2], date[1],
				units="hours"))
			y = data.frame(lon=y[,1], lat=y[,2])
			},
		DCRW3S = {
			b = apply(psamples$b, 1, mean)
			b.5 = apply(psamples$b, 1, median)
			summary = data.frame(id=as.character(id), date =
				as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
				lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
				lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,], b, b.5)
			model = model
			mcmc.settings = list(burnin = adapt, posterior.samples = samples, thinning = thin, n.chains = chains)
			step = with(summary, difftime(date[2], date[1],
				units="hours"))
			y = data.frame(lon=y[,1], lat=y[,2])
			})
	out = list(summary=summary, mcmc=psamples, model=model, mcmc.settings=mcmc.settings,
		timestep=step, N=RegN, data=y)
	out
    }
    lapply(loc.list, ssm1)
}
