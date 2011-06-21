\name{CDMN0P.Lik}
\alias{CDMN0P.Lik}
\title{
Likelihood Evaluation of Pure Depletion Catch Dynamic Model
}
\description{
This function evaluates a modified profile likelihood approximation to the likelihood
of the catch data as a function parameter values when there is no in-season perturbation.
}
\usage{
CDMN0P.Lik(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
A numeric vector of length 5 with initial parameter values in the log scale
}
  \item{dates}{
An integer vector of length 2 wih initial and final time steps.
}
  \item{obscat}{
A numeric vector with observed catch in numbers by time step..
}
  \item{obseff}{
A numeric vector with observed nominal effort by time step..
}
  \item{M.fixed}{
Logical, shall natural be estimated (M.fixed=FALSE) or remain fixed at a given value?
}
  \item{M}{
If M.fixed is TRUE, this number sets the value of the natural mortality rate.
}
  \item{distr}{
Character. 'normal' or 'lognormal' model for the observed random variable, the catch.
}
}
\details{
This function defines a pure depletion model to be fit, with no in-season perturbations.
Probably good as a basic -null- model.

If the catch shows an increasing trend over the season this model may yield convergence problems.
In that case the basic model may be a model with one or two perturbations.
}
\value{
A scalar, the modified profile likelihood evaluated at given parameter values.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Parameter values are to be in log scale to avoid scaling issues during the
numerical optimization. Back-transformed MLEs and their standard errors are computed
by the wrapper function catdyn().
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 8e-6
N0.ini.B.1990.S1        <- 2.6
k.ini.B.1990.S1         <- 3e-4
alpha.ini.B.1990.S1     <- 1.32
beta.ini.B.1990.S1      <- 0.45
pars.ini.B.1990.S1.0P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
dates.B.1990.S1.0P      <- c(head(SeasonData.1990.S1.B$period,1),
                             tail(SeasonData.1990.S1.B$period,2))
#Function evaluation at initial parameter values
CDMN0P.Lik(par=pars.ini.B.1990.S1.0P,
           dates=dates.B.1990.S1.0P,
           obscat=SeasonData.1990.S1.B$obscat,
           obseff=SeasonData.1990.S1.B$obseff2,
           M.fixed=FALSE,
           distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }

