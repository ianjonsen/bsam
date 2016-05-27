`fitSSM` <-
function (indata, model="DCRW", tstep=1, adapt=20000, samples=10000, thin=10, chains=2)
{
	if(!model %in% c('DCRW', 'DCRWS', 'DCRW3S', 'hDCRW', 'hDCRWS')) stop("Model not implemented")
    sf = system.file(package = "bsam")
    if (.Platform[1] == "unix") {
 		model.file = paste(sf, "/jags/", model, ".txt", sep = "")
    	}  	
    else if (.Platform[1] == "windows") {
        model.file = paste(sf, "\\jags\\", model, ".txt", 
            sep = "")
    	}

	options(warn=-1)	    	
    seed = sample(1:1e+05, 1)
    st = proc.time()
      
	data = dat4jags(indata, tstep=tstep, tod=TRUE)	
	switch(model, 
		DCRW = {
			fit = ssm(data, model = model, model.file = model.file, adapt = adapt, samples = samples, thin = thin, chains = chains)
		},
		DCRWS = {
			fit = ssm(data, model = model, model.file = model.file, adapt = adapt, samples = samples, thin = thin, chains = chains)
		},
		DCRW3S = {
			fit = ssm(data, model = model, model.file = model.file, adapt = adapt, samples = samples, thin = thin, chains = chains)
		},	
		hDCRW = {	
			fit = hssm(data, model = model, model.file = model.file, adapt = adapt, samples = samples, thin = thin, chains = chains)
		},
		hDCRWS = {	
			fit = hssm(data, model = model, model.file = model.file, adapt = adapt, samples = samples, thin = thin, chains = chains)
		})		
	cat("Elapsed time: ", round((proc.time() - st)[3]/60,2), "min \n")	
	options(warn=0)
	
	fit
}
