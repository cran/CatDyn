\name{CDMN3P.Lik}
\alias{CDMN3P.Lik}
\title{
Likelihood Evaluation of 3-Perturbations Catch Dynamic Model
}
\description{
This function evaluates a modified profile likelihood approximation to the likelihood
of the catch data as a function parameter values when there are three in-season perturbations.
}
\usage{
CDMN3P.Lik(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
Numeric vector of length 8 with initial parameter values in log scale.
}
  \item{dates}{
An integer vector of length 5 with initial time step, perturbation time steps, and final time step.
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
This function defines a 3 perturbations model to be fit.
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
M                       <- 2e-3
N0.ini.B.1990.S1        <- 4.3
P1.ini.B.1990.S1        <- 0.9
P2.ini.B.1990.S1        <- 0.15
P3.ini.B.1990.S1        <- 0.3
k.ini.B.1990.S1         <- 8e-6
alpha.ini.B.1990.S1     <- 1.7
beta.ini.B.1990.S1      <- 2.4
pars.ini.B.1990.S1.3P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(P1.ini.B.1990.S1),
                             log(P2.ini.B.1990.S1),
                             log(P3.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
P1.B.1990.S1.3P         <- 71
P2.B.1990.S1.3P         <- 86
P3.B.1990.S1.3P         <- 135
dates.B.1990.S1.3P      <- c(head(SeasonData.1990.S1.B$period,1),
                             P1.B.1990.S1.3P,
                             P2.B.1990.S1.3P,
                             P3.B.1990.S1.3P,
                             tail(SeasonData.1990.S1.B$period,1))
#Function evaluation at initial parameter values
CDMN3P.Lik(par=pars.ini.B.1990.S1.3P,
           dates=dates.B.1990.S1.3P,
           obscat=SeasonData.1990.S1.B$obscat,
           obseff=SeasonData.1990.S1.B$obseff2,
           M.fixed=FALSE,
           distr='lognormal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }