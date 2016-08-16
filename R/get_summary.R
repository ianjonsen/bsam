#' Extract summary output and optionally export as a .csv file.
#' 
#' Takes a fitted \code{fit_ssm} object and extracts the summary data.frame, which includes
#' the animal ids, POSIXct date/time (at increments specified by \code{tstep} in the \code{fit_ssm} call), 
#' posterior mean longitude and latitude, and the 2.5, 50, and 97.5 %-iles of the posterior MCMC samples for 
#' longitude and latitude. For the \code{DCRWS} and \code{hDCRWS} models, the posterior mean and median behavioural
#' states corresponding to each estimated location are also provided.
#' 
#' @param x an output object from \code{fit_ssm}. If not an error will be returned.
#' @param file a character string naming a file. " " indicates output to the console (default)
#' @return a summary data.frame printed either to the console (default) or written as .csv to a
#' specified file.
#' 
#' @examples
#' \dontrun{
#' data(ellie)
#' fit <- fit_ssm(ellie, model = "DCRW", tstep = 1, adapt = 2000, samples = 1000, 
#'                 thin = 2, span = 0.1)
#'                 
#' ## print to console                 
#' get_summary(fit)
#' 
#' ## export to .csv file
#' get_summary(fit, file = "ellie_dcrw.csv")
#' 
#' ## plot track of second animal (posterior mean longitude and latitude)
#' p <- ggplot(data = subset(get_summary(fit), id == unique(id)[2])) + geom_point(aes(lat, lon), col = "firebrick", size = 2)
#' p
#' 
#' ## overlay posterior median longitude and latitude
#' p + geom_point(aes(lat.5, lon.5), col = "dodgerblue", size = 1)
#' }
#' @importFrom utils write.csv
#' @export 

get_summary <- function(x, file = " ") 
{
  if(class(x) != "ssm" && class(x) != "hssm") stop("Input is not a fit_ssm object")

  else if(class(x) == "ssm") {
    summ <- do.call(rbind, lapply(x, function(z) z$summary))
  }
  
  else if(class(x) == "hssm") {
    summ <- x$summary
  }

  if(file == " ") return(summ)
  
  else {
    write.csv(summ, file = file, row.names = FALSE)
  }
}