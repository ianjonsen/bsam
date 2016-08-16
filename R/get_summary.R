#' Extract summary output as a \code{dplyr::data_frame} and optionally export as a .csv file.
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
#' get_summary(fit, file = "dcrw.csv")
#' }
#' @importFrom dplyr as_data_frame
#' @export 
#' 
get_summary <- function(x, file = " ") 
{
  summ <- dplyr::as_data_frame(do.call(rbind, lapply(x, function(z) z$summary)))

  if(file == " ") return(summ)
  
  else {
    write.csv(summ, file = file, row.names = FALSE)
  }
}