\name{CatDynPar}
\alias{CatDynPar}
\title{
Table of Parameter Estimates and Coefficients of Variation
}
\description{
It takes parameter estimates from output by CatDynFit() from a single model fit
and creates a neat table of estimates, dates in ISO 8601 format for each abundance
input/output, and coefficient of variation of estimates.
}
\usage{
CatDynPar(x, method, partial = TRUE)
}
\arguments{
  \item{x}{
A single model object of class catdyn output by CatDynFit().
}
  \item{method}{
The numerical method used to minimize the loglikelihod of x.
}
  \item{partial}{
Logical, FALSE for transit models, TRUE (the default) for all other model versions.
}
}
\details{
Arguments x and method must both be of length = 1.
}
\value{
A data.frame with the following columns: parameter names (all estimated parameters,
including the dispersion parameter for likelihood options other than Poisson, adjusted
profile normal, and adjusted profile lognormal), timing (the data in ISO 8601 format
of each input/output event), the maximum likelihood estimates, and the percent
coefficient of variation. In two fleet models the last three columns are duplicated,
with the fleet name added to the column names.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
In model fits with a monthly time step, the timing column includes year and month.
In model fits with a weekly time step, the timing includes the full date with initial
and final date of that week.

In model fits with a daily time step, the timing includes the full data.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models }