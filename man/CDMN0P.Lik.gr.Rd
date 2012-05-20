\name{CDMN0P.Lik.gr}
\alias{CDMN0P.Lik.gr}
\title{
Calculate Analytical Gradients for the Pure Depletion Model
}
\description{
Given parameter values, the function calculates analytical gradients for each of
five dimensions of the likelihood function of the pure depletion model from
formulas obtained using the Maxima Computer Algebra System 5.24.00.
}
\usage{
CDMN0P.Lik.gr(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
Numeric vector of length 5 with initial parameter values in log scale.
}
  \item{dates}{
An integer vector of length 2 wih initial and final time steps.
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
Character, either 'normal' or 'lognormal'.
}
}
\details{
The likelihood approximation is the modified profile likelihood function.
}
\value{
A numeric vector of gradients of length 5 for given parameter values.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Parameter values are to be in log scale for consistency with CDMN0P.Lik.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 0.001
N0.ini.B.1990.S1        <- 5.3
k.ini.B.1990.S1         <- 0.00027
alpha.ini.B.1990.S1     <- 1.32
beta.ini.B.1990.S1      <- 0.45
pars.ini.B.1990.S1.0P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
dates.B.1990.S1.0P      <- c(head(SeasonData.1990.S1.B$period,1),
                             tail(SeasonData.1990.S1.B$period,1))
#Gradients at the initial parameter values
CDMN0P.Lik.gr(par=pars.ini.B.1990.S1.0P,
              dates=dates.B.1990.S1.0P,
              obscat=SeasonData.1990.S1.B$obscat,
              obseff=SeasonData.1990.S1.B$obseff2,
              M.fixed=FALSE,
              distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }

