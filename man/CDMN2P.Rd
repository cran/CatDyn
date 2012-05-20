\name{CDMN2P}
\alias{CDMN2P}
\title{
Implementing a 2-Perturbations Catch Dynamics Model to Examine Parameter Values
}
\description{
This function will do exactly the same as its homonym CDMN2P.Lik, except that instead
of producing an evaluated approximation to the likelihood, it outputs a data.frame
with several results of interest.
}
\usage{
CDMN2P(par, dates, obscat, obseff, obsmbm, M.fixed, M, distr)
}
\arguments{
  \item{par}{
A numeric vector of length 7 with initial parameter values in the log scale.
}
  \item{dates}{
An integer vector of length 4 with initial time step, perturbation time steps, and final time step.
}
  \item{obscat}{
A numeric vector with observed catch in numbers by time step.
}
  \item{obseff}{
A numeric vector with observed nominal effort by time step.
}
  \item{obsmbm}{
A numeric vector with observed mean body mass by time step.
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
Parameter values are to be in log scale.
}
\value{
A data.frame, an object of class CatDynMod, with the following columns:
the period, the sequence of time steps defining the season,
the observed effort,
the observed catch, in billions,
the model predicted catch, in billions,
the residuals,
the biomass in tonnes.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Thus is used to feed a standard plot for exploratory analysis.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 0.01
N0.ini.B.1990.S1        <- 5.3
P1.ini.B.1990.S1        <- 1.5
P2.ini.B.1990.S1        <- 0.5
k.ini.B.1990.S1         <- 0.00027
alpha.ini.B.1990.S1     <- 1.32
beta.ini.B.1990.S1      <- 0.45
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
#Results at initial parameter values
CDMN2P(par=pars.ini.B.1990.S1.2P,
       dates=dates.B.1990.S1.2P,
       obscat=SeasonData.1990.S1.B$obscat,
       obseff=SeasonData.1990.S1.B$obseff2,
       obsmbm=SeasonData.1990.S1.B$obsmbm,
       M.fixed=FALSE,
       distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
