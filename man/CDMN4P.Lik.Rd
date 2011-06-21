\name{CDMN4P.Lik}
\alias{CDMN4P.Lik}
\title{
Likelihood Evaluation of 4-Perturbations Catch Dynamic Model
}
\description{
This function evaluates a modified profile likelihood approximation to the likelihood
of the catch data as a function parameter values when there are four in-season perturbations.
}
\usage{
CDMN4P.Lik(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
Numeric vector of length 9 with initial parameter values in log scale.
}
  \item{dates}{
An integer vector of length 6 with initial time step, perturbation time steps, and final time step.
}
  \item{obscat}{
A numeric vector with observed catch in numbers by time step.
}
  \item{obseff}{
A numeric vector with observed nominal effort by time step.
}
  \item{M.fixed}{
Logical, shall natural be estimated (M.fixed=FALSE) or remain fixed at a given value?
}
  \item{M}{
If M.fixed is TRUE, this number sets the value of the natural mortality rate.
}
  \item{distr}{
Character, either 'normal' or 'lognormal'.
}
}
\details{
This function defines a 4 perturbations model to be fit.
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
M                       <- 1.3e-2
N01.ini.B.1990.S1       <- 37.1
P1.ini.B.1990.S1        <- 2.0
P2.ini.B.1990.S1        <- 10.6
P3.ini.B.1990.S1        <- 0.1
P4.ini.B.1990.S1        <- 4.7
k.ini.B.1990.S1         <- 8e-6
alpha.ini.B.1990.S1     <- 1.4
beta.ini.B.1990.S1      <- 0.95
pars.ini.B.1990.S1.4P   <- c(log(M),
                             log(N01.ini.B.1990.S1),
                             log(P1.ini.B.1990.S1),
                             log(P2.ini.B.1990.S1),
                             log(P3.ini.B.1990.S1),
                             log(P4.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
P1.B.1990.S1.4P         <- 56
P2.B.1990.S1.4P         <- 71
P3.B.1990.S1.4P         <- 86
P4.B.1990.S1.4P         <- 135
dates.B.1990.S1.4P      <- c(head(SeasonData.1990.S1.B$period,1),
                                  P1.B.1990.S1.4P,
                                  P2.B.1990.S1.4P,
                                  P3.B.1990.S1.4P,
                                  P4.B.1990.S1.4P,
                                  tail(SeasonData.1990.S1.B$period,1))
sealen.B.1990.S1        <- dates.B.1990.S1.4P[2]-dates.B.1990.S1.4P[1]+1
#Function evaluation at initial parameter values
CDMN4P.Lik(par=pars.ini.B.1990.S1.4P,
           dates=dates.B.1990.S1.4P,
           obscat=SeasonData.1990.S1.B$obscat,
           obseff=SeasonData.1990.S1.B$obseff2,
           M.fixed=FALSE,
           distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }