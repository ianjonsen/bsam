\name{dat4jags}
\alias{dat4jags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepares and writes input data for DCRW/DCRWS/hDCRWS models}

\description{
  Takes an R data.frame of Argos tracking data restructures it for the models. Intended for internal use, called by fitSSM.
}
\usage{
dat4jags(locfile, tstep = 1, tod = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locfile}{An R data.frame of Argos tracking data. See fitSSM (indata) for details on structure.}
  \item{tstep}{The time step to be assumed for the state-space model, specified in days.}
  \item{tod}{Logical. Specifies if absolute time of day of observations is relevant.}
  \item{\dots}{Other arguments may be passed.}
}
\value{
  Returns a list to be used by ssm/hssm.

}

\seealso{Function to be called by \code{\link{fitSSM}}.}
