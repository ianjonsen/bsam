#' @param loc.list 
#'
#' @param model 
#' @param adapt 
#' @param samples 
#' @param thin 
#' @param chains 
#' @param ... 
#'
#' @export
#' @importFrom rjags jags.model
`hssm` =
function (loc.list, model = "hDCRWS", adapt, samples, thin, chains, ...)
{
    if (!model %in% c("hDCRW","hDCRWS"))
        stop("model not implemented.")

    nt = length(loc.list)
    y = NULL
    itau2 = NULL
    nu = NULL
    j = NULL
    RegN = NULL
    id = NULL
    first.date = NULL
    tstep = NULL
    len.y = sapply(loc.list, function(x) nrow(x$y))
    ids = sapply(loc.list, function(x) x$id)
    for(i in 1:nt){
    	y = rbind(y, loc.list[[i]]$y)
    	itau2 = rbind(itau2, loc.list[[i]]$itau2)
    	nu = rbind(nu, loc.list[[i]]$nu)
    	j = c(j, loc.list[[i]]$j)
    	RegN = c(RegN, loc.list[[i]]$RegN)
    	id = c(id, loc.list[[i]]$id)
    	first.date = c(first.date, loc.list[[i]]$first.date)
    	tstep = c(tstep, loc.list[[i]]$tstep)
    	}
    tmp = list(NULL)
    for(i in 1:nt){
		tmp[[i]] = loc.list[[i]]$idx
    	}
 	for(i in 2:nt){
 		tmp[[i]] = max(tmp[[i-1]]) + tmp[[i]][-1]-1
 		}

 	extract = function(k) nrow(k$y)

 	idx = unlist(tmp)
	Xidx = cumsum(c(1, RegN-1))
	# Yidx can not be used in JAGS implementation but leave it here for now
	Yidx = cumsum(c(1, unlist(lapply(loc.list, extract))))

    x = rbind(y[idx[-length(idx)], ], y[max(idx) - 1, ])
    row.na = which(is.na(x[, 1]))
    x[row.na, 1] = approx(seq(nrow(x)), x[, 1], xout = row.na,
            rule = 2)$y
    row.na = which(is.na(x[, 2]))
    x[row.na, 2] = approx(seq(nrow(x)), x[, 2], xout = row.na,
            rule = 2)$y
    x = x[-nrow(x),]
    row.na = which(is.na(y[, 1]) | is.na(y[, 2]))
    if (length(row.na) == 0) {
        y.init = y
        y.init[!is.na(y.init)] = NA
    	}
    else {
       y.init = y
       y.init[!is.na(y.init)] = 9999
       yinit.xna = approx(seq(nrow(y)), y[, 1], xout = row.na,
                rule = 2)$y
       yinit.yna = approx(seq(nrow(y)), y[, 2], xout = row.na,
                rule = 2)$y
       y.init[row.na, ] = cbind(yinit.xna, yinit.yna)
       y.init[y.init == 9999] = NA
       }
    tstep.sec = tstep[1] * 86400 #only need first element; all tsteps the same
    lidx.fn = function(k) length(k$idx) - 1
    lidx = sapply(loc.list, lidx.fn)
    steplims = unlist(sapply(1:nt, function(i){
    				seq(first.date[i], by = tstep.sec,
            		length = lidx[i])
            		}))

	jags.data = list(y = y, itau2 = itau2, nu = nu, j = j, idx = idx, Xidx=Xidx,
		Yidx=Yidx, N=nt)
    iSigma = matrix(c(1, 0, 0, 1), 2, 2)

	switch(model,
        hDCRW = {
            jags.inits = list(list(iSigma = iSigma, gamma = 0.5,
             theta = 0, x = x, logpsi = runif(nt,-1,1), y = y.init),
             list(iSigma = iSigma, gamma = 0.4, theta = 0.1, x = x,
             logpsi = runif(nt,-1,1), y = y.init), list(iSigma = iSigma, gamma = 0.45,
             theta = -0.1, x = x, logpsi = runif(nt,-1,1), y = y.init),
             list(iSigma = iSigma, gamma = 0.3, theta = -0.2, x = x,
             logpsi = runif(nt,-1,1), y = y.init))
            jags.params = c("Sigma", "x", "theta", "gamma", "psi")
        	},
		hDCRWS = {
		    jags.inits = list(list(iSigma = iSigma, gamma = c(0.8, NA), dev = 0.6,
    		 tmp = c(0.45, 0.55), alpha = c(0.55, 0.45), lambda = c(0.45, NA),
    		 b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(nt,-1,1), y = y.init),
    		 list(iSigma = iSigma, gamma = c(0.85, NA), dev = 0.65, tmp = c(0.6, 0.4),
    		 alpha = c(0.4, 0.6), lambda = c(0.5, NA), b = rbinom(nrow(x), 1, 0.5) + 1,
    		 x = x, logpsi = runif(nt,-1,1), y = y.init), list(iSigma = iSigma,
    		 gamma = c(0.75, NA), dev = 0.7, tmp = c(0.5, 0.6), alpha = c(0.5, 0.5),
    		 lambda = c(0.55, NA), b = rbinom(nrow(x), 1, 0.5) + 1, x = x,
    		 logpsi = runif(nt,-1,1), y = y.init), list(iSigma = iSigma, gamma = c(0.6, NA),
    		 dev = 0.5, tmp = c(0.6, 0.5), alpha = c(0.65, 0.45), lambda = c(0.4, NA),
    		 b = rbinom(nrow(x), 1, 0.5) + 1, x = x, logpsi = runif(nt,-1,1), y = y.init))
    		jags.params = c("Sigma", "x", "theta", "gamma", "alpha", "b", "psi")
	    	})

	if(chains==1) jags.inits = jags.inits[[1]]
	if(chains==2) jags.inits = list(jags.inits[[1]],jags.inits[[2]])
	if(chains==3) jags.inits = list(jags.inits[[1]],jags.inits[[2]],jags.inits[[3]])

	model.file = paste(system.file('jags',package='bsam'), "/", model, ".txt", sep="")

	burn = jags.model(model.file, jags.data, jags.inits, n.chains=chains, n.adapt=adapt/2)
	update(burn, n.iter=adapt/2)
	psamples = jags.samples(burn, jags.params, n.iter=samples, thin=thin)

	lon = apply(psamples$x[,1,,],1,mean)
	lat = apply(psamples$x[,2,,],1,mean)
	lon.q = apply(psamples$x[,1,,],1, quantile, c(0.025, 0.5, 0.975))
	lat.q = apply(psamples$x[,2,,],1, quantile, c(0.025, 0.5, 0.975))
	if(model=="hDCRWS"){
		b = apply(psamples$b,1,mean)
		b.5 = apply(psamples$b,1,median)
		}
	switch(model,
		hDCRW = {
			summary = data.frame(id=rep(as.character(ids), RegN-1), date =
				as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
				lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
				lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,])
			},
		hDCRWS = {
			summary = data.frame(id=rep(as.character(ids), RegN-1), date =
				as.POSIXct(as.numeric(steplims), origin="1970-01-01 00:00:00", tz="GMT"),
				lon, lat, lon.025=lon.q[1,], lon.5=lon.q[2,], lon.975=lon.q[3,],
				lat.025=lat.q[1,], lat.5=lat.q[2,], lat.975=lat.q[3,], b, b.5)
			})

	model = model
	mcmc.settings = list(burnin = adapt, posterior.samples = samples, thinning = thin,
		n.chains = chains)
	step = with(summary, difftime(date[2], date[1], units="hours"))
	y = data.frame(id=rep(ids,len.y),lon=y[,1], lat=y[,2])
	N = RegN - 1
	numtracks = length(RegN)
	out = list(summary=summary, mcmc=psamples, model=model, mcmc.settings=mcmc.settings,
		timestep=step, N=N, data=y, numtracks=numtracks)

	out
}
