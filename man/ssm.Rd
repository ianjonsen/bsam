\name{ssm}
\alias{ssm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits state-space model to Argos data}
\description{
  Takes output from dat4jags, sets up initial values, calls JAGS, and aggregates results. Intended for internal use, called by fitSSM.}
\usage{
ssm(loc.list, model = "DCRW", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc.list}{Data from dat4jags to be passed to JAGS}
  \item{model}{The state-space model to be fit. Currently, DCRW and DCRWS are the only model options.}
  \item{\dots}{Other arguments can be passed.}
}

\value{Returns a list of McMC samples from marginal posteriors and a summary data.frame of mean and median position estimates.}
\seealso{Function to be called by \code{\link{fitSSM}}.}
