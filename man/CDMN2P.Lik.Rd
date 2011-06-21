\name{CDMN2P.Lik}
\alias{CDMN2P.Lik}
\title{
Likelihood Evaluation of 2-Perturbations Catch Dynamic Model
}
\description{
This function evaluates a modified profile likelihood approximation to the likelihood
of the catch data as a function parameter values when there are two in-season perturbations.
}
\usage{
CDMN2P.Lik(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
A numeric vector of length 7 initial parameter values in the log scale
}
  \item{dates}{
An integer vector of length 4 with initial time step, perturbation time steps, and final time step.
}
  \item{obscat}{
A numeric vector with observed catch in numbers by time step.
}
  \item{obseff}{
A numeric vector with observed daily nominal effort by time step.
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
This function defines a 2 perturbations model to be fit.
}
\value{
A scalar, the modified profile likelihood evaluated at given parameter values.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Parameter values are to be in log scale to avoid scaling issues during the
numerical minimization. Back-transformed MLEs and their standard errors are computed
by the wrapper function catdyn().
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 4e-3
N0.ini.B.1990.S1        <- 8.5
P1.ini.B.1990.S1        <- 1.4
P2.ini.B.1990.S1        <- 0.5
k.ini.B.1990.S1         <- 4e-6
alpha.ini.B.1990.S1     <- 1.4
beta.ini.B.1990.S1      <- 3.1
pars.ini.B.1990.S1.2P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(P1.ini.B.1990.S1),
                             log(P2.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
P1.B.1990.S1.2P         <- 71
P2.B.1990.S1.2P         <- 135
dates.B.1990.S1.2P      <- c(head(SeasonData.1990.S1.B$period,1),
                             P1.B.1990.S1.2P,
                             P2.B.1990.S1.2P,
                             tail(SeasonData.1990.S1.B$period,1))
#Function evaluation at initial parameter values
CDMN2P.Lik(par=pars.ini.B.1990.S1.2P,
           dates=dates.B.1990.S1.2P,
           obscat=SeasonData.1990.S1.B$obscat,
           obseff=SeasonData.1990.S1.B$obseff2,
           M.fixed=FALSE,
           distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }

