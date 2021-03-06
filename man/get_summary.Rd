% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_summary.R
\name{get_summary}
\alias{get_summary}
\title{Extract summary output and optionally export as a .csv file.}
\usage{
get_summary(x, file = " ")
}
\arguments{
\item{x}{an output object from \code{fit_ssm}. If not an error will be returned.}

\item{file}{a character string naming a file. " " indicates output to the console (default)}
}
\value{
a summary data.frame printed either to the console (default) or written as .csv to a
specified file.
}
\description{
Takes a fitted \code{fit_ssm} object and extracts the summary data.frame, which includes
the animal ids, POSIXct date/time (at increments specified by \code{tstep} in the \code{fit_ssm} call), 
posterior mean longitude and latitude, and the 2.5, 50, and 97.5 %-iles of the posterior MCMC samples for 
longitude and latitude. For the \code{DCRWS} and \code{hDCRWS} models, the posterior mean and median behavioural
states corresponding to each estimated location are also provided.
}
