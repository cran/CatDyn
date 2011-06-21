\name{CDMN2P.Lik.gr}
\alias{CDMN2P.Lik.gr}

\title{
Calculate Analytical Gradients for the 2-Perturbations Model (Beta Testing)
}
\description{
Given parameter values, the function calculate analytical gradients for each of
the seven dimensions of the objective function of the 2-Perturbations model from
formulas obtained using the Maxima Computer Algebra System 5.24.00.
}
\usage{
CDMN2P.Lik.gr(par, dates, obscat, obseff, M.fixed, M, distr)
}

\arguments{
  \item{par}{
Numeric vector of length 7 with initial parameter values in log scale.
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
The likelihood approximation is the modified profile likelihood function.
}
\value{
A numeric vector of gradients of length 7 for given parameter values.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Parameter values are to be in log scale for consistency with CDMN2P.Lik.
Checked against numerical gradients obtained with optimx, method 'spg'.
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
#Gradients at the initial parameter values
CDMN2P.Lik.gr(par=pars.ini.B.1990.S1.2P,
              dates=dates.B.1990.S1.2P,
              obscat=SeasonData.1990.S1.B$obscat,
              obseff=SeasonData.1990.S1.B$obseff2,
              M.fixed=FALSE,
              distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
