\name{CDMN1P.Lik}
\alias{CDMN1P.Lik}

\title{
Likelihood Evaluation of 1-Perturbation Catch Dynamics Model
}
\description{
This function evaluates a modified profile likelihood approximation to the likelihood
of the catch data as a function parameter values when there is one in-season perturbation.
}
\usage{
CDMN1P.Lik(par, dates, obscat, obseff, M.fixed, M, distr)
}

\arguments{
  \item{par}{
A numeric vector of length 6 with initial parameter values in the log scale
}
  \item{dates}{
An integer vector of length 3 with initial time step, perturbation time step, and final time step.
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
A positive number or NULL; if M.fixed is TRUE, this number sets the value of the natural mortality rate.
}
  \item{distr}{
Character, 'normal' or 'lognormal' model for the observed random variable, the catch.
}
}
\details{
This function defines a 1 perturbation model to be fit.
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
M                       <- 5e-5
N0.ini.B.1990.S1        <- 5.1
P1.ini.B.1990.S1        <- 1.5
k.ini.B.1990.S1         <- 7e-7
alpha.ini.B.1990.S1     <- 1.7
beta.ini.B.1990.S1      <- 3.1
pars.ini.B.1990.S1.1P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(P1.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
P1.B.1990.S1.1P         <- 71
dates.B.1990.S1.1P      <- c(head(SeasonData.1990.S1.B$period,1),
                             P1.B.1990.S1.1P,
                             tail(SeasonData.1990.S1.B$period,1))
#Function evaluation at initial parameter values
CDMN1P.Lik(par=pars.ini.B.1990.S1.1P,
           dates=dates.B.1990.S1.1P,
           obscat=SeasonData.1990.S1.B$obscat,
           obseff=SeasonData.1990.S1.B$obseff2,
           M.fixed=FALSE,
           distr='lognormal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }

