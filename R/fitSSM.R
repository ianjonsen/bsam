#' Fit Bayesian state-space models to animal movement data
#' 
#' Fits state-space models to Argos satellite tracking data. User can choose
#' between a first difference correlated random walk (DCRW) model, a switching
#' model (DCRWS), and a hierarchical swtiching model (hDCRWS) for estimating
#' behavioural states.
#' 
#' The models are fit using JAGS 3.1.0 (Just Another Gibbs Sampler, created and
#' maintained by Martyn Plummer; http://martynplummer.wordpress.com/;
#' http://mcmc-jags.sourceforge.net). fitSSM is a wrapper function that first
#' calls dat4jags(), which prepares the input data, then calls ssm()/hssm(),
#' which fits the specified state-space model to the data, returning a list of
#' results.
#' 
#' @param indata An R data.frame containing the following columns, "id","date",
#' "lc", "lon", "lat". "id" is a unique identifier for the tracking dataset.
#' "date" is the GMT date-time of each observation with the following format
#' "2001-11-13 07:59:59". "lc" is the Argos location quality class of each
#' observation, values in ascending order of quality are "B", "A", "0", "1",
#' "2", "3". "lon" is the observed longitude in decimal degrees. "lat" is the
#' observed latitude in decimal degress.
#' @param model name of state-space model to be fit to data. Currently, this
#' can be one of "DCRW", "DCRWS", or "hDCRWS".
#' @param tstep time step as fraction of a day, default is 1 (24 hours).
#' @param adapt number of samples during the adaptation and update (burn-in)
#' phase, adaptation and updates are fixed at adapt/2
#' @param samples number of samples to generate after convergence is assumed
#' @param thin amount of thining of to be applied to minimize sample
#' autocorrelation.
#' @return For DCRW and DCRWS models, a list is returned with each outer list
#' elements corresponding to each unique individual id in the input data.
#' Within these outer elements are a "summary" data.frame of posterior mean and
#' median state estimates (locations or locations and behavioural states), the
#' name of the "model" fit, the "timestep" used, the input location "data", the
#' number of location state estimates ("N"), and the full set of "mcmc"
#' samples. For the hDCRWS model, a list is returned where results, etc are
#' combined amongst the individuals.
#' @author Ian Jonsen
#' @references Jonsen ID, Myers RA, Mills Flemming J (2003) Meta-analysis of
#' animal movement using state-space models. Ecology 84:3055-3063
#' 
#' Jonsen ID, Mills Flemming J, Myers RA (2005) Robust state-space modeling of
#' animal movement data. Ecology 86:2874-2880
#' 
#' Block et al. (2011) Tracking apex marine predator movements in a dynamic
#' ocean. Nature 475:86-90
#' 
#' Jonsen et al. (2013) State-space models for biologgers: a methodological
#' road map. Deep Sea Research II DOI: 10.1016/j.dsr2.2012.07.008
#' 
#' @examples
#' # Fit DCRW model for state filtering and regularization
#' data(lbt)
#' #fit <- fitSSM(lbt, model="DCRW", tstep=1, adapt=40000, samples=20000, thin=20)
#' #plotSSM(fit, save.to.pdf=FALSE)
#' #diagSSM(fit, save.to.pdf=FALSE)
#' 
#' # Fit DCRWS model for state filtering, regularization and behavioural state estimation
#' # Not run
#' # data(lbt)
#' # fit <- fitSSM(lbt, model="DCRWS", tstep=0.5, adapt=40000, samples=20000, thin=20)
#' # plotSSM(fit, save.to.pdf=FALSE)
#' # diagSSM(fit, save.to.pdf=FALSE)
#' 
#' # fit hDCRWS model to >1 tracks simultaneously
#' # this can provide better parameter and behavioural state estimation 
#' # by borrowing strength across multiple track datasets
#' 
#' # Not run
#' # data(lbt)
#' # tmp <- lbt
#' # tmp$id <- 15395
#' # lbt2 <- rbind(lbt,tmp)
#' # This will take some time to complete
#' # fit <- fitSSM(lbt2, model="hDCRWS", tstep=0.5, adapt=40000, samples=20000, thin=20)
#' # plotSSM(fit, save.to.pdf=FALSE)
#' # diagSSM(fit, save.to.pdf=FALSE)
#' 
#' @export 
`fitSSM` <-
function (d, model="DCRW", tstep=1, adapt=40000, samples=20000, thin=20, span=0.2)
{
	if(!model %in% c('DCRW', 'DCRWS', 'hDCRW', 'hDCRWS')) stop("Model not implemented")
  model.file <- file.path(system.file(package = "bsam"), "jags", paste(model, ".txt", sep=""))
    
	options(warn=-1)	    	
  seed <- sample(1:1e+05, 1)
  st <- proc.time()
      
	dd <- dat4jags(d, tstep = tstep, tpar=tpar())	
	if(model %in% c("DCRW","DCRWS")) {
	  fit <- ssm(dd, model = model, adapt = adapt, samples = samples, thin = thin, 
	             chains = 2, span = span)
	}
	else {
	  fit <- hssm(dd, model = model, adapt = adapt, samples = samples, thin = thin, 
	              chains = 2, span = span)
	}
	
	cat("Elapsed time: ", round((proc.time() - st)[3]/60,2), "min \n")	
	options(warn=0)
	
	fit
}
