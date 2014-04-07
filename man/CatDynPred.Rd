\name{CatDynPred}
\alias{CatDynPred}
\title{
Predictions for a Catch Dynamic Model
}
\description{
Calculates the predicted catch, residuals, by fleet, and population biomass of a 
fish stock using parameters estimated for a Catch Dynamic Model.
}
\usage{
CatDynPred(x, method)
}
\arguments{
  \item{x}{
An object of class catdyn.
}
  \item{method}{
Character. The particular numerical method from which estimates have to be drawn.
}
}
\details{
This function is very similar to catdynexp() but instead of using arbitrary
parameter values given by the user, it takes the maximum likelihood estimates
produced by the optimizer.
}
\value{
A list of length 2.
\item{Properties }{A list of length 3. 
'Units' is a dataframe with the units of time step, 
catch, body mass, and the numbers multiplier. 'Fleets' is a dataframe with the fleets
names and the units of nominal effort for each fleet. 'Dates' is a dataframe with start 
and end dates of the fishing season or period in the ISO 8601 format}
\item{Model }{A list of length 5.
'Type' is the perturbation type of model. 'Dates' is the timing of perturbations,
'Distr' is the chosen distribution for the observed catch data, 'Parameters' are the 
maximum likelihood estimates from CatDynFit(), and 'Results', is a dataframe with the time step, 
and for each fleet, the observed effort, observed catch, predicted catch, and the 
residuals, plus two more columns, one with the predicted population abundance, 
and a final one with the predicted population biomass }
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models  }

