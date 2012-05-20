\name{CDMN1P.Lik.gr}
\alias{CDMN1P.Lik.gr}
\title{
Calculate Analytical Gradients for the 1-Perturbation Model
}
\description{
Given parameter values, the function calculates analytical gradients for each of
the six dimensions of the likelihood function of the 1-Perturbation model from
formulas obtained using the Maxima Computer Algebra System 5.24.00.
}
\usage{
CDMN1P.Lik.gr(par, dates, obscat, obseff, M.fixed, M, distr)
}
\arguments{
  \item{par}{
Numeric vector of length 6 with initial parameter values in log scale.
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
Character, either 'normal' or 'lognormal' for the distribution of the response, catch.
}
}
\details{
The likelihood approximation is the modified profile likelihood function.
}
\value{
A numeric vector of gradients of length 6 for given parameter values.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Parameter values are to be in log scale for consistency with CDMN1P.Lik.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 0.01
N0.ini.B.1990.S1        <- 5.3
P1.ini.B.1990.S1        <- 0.5
k.ini.B.1990.S1         <- 0.00027
alpha.ini.B.1990.S1     <- 1.32
beta.ini.B.1990.S1      <- 0.45
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
#Gradients at the initial parameter values
CDMN1P.Lik.gr(par=pars.ini.B.1990.S1.1P,
              dates=dates.B.1990.S1.1P,
              obscat=SeasonData.1990.S1.B$obscat,
              obseff=SeasonData.1990.S1.B$obseff2,
              M.fixed=FALSE,
              distr='normal')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }

