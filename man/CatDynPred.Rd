\name{CatDynPred}
\alias{CatDynPred}
\title{
Predictions from a Generalized Depletion (Catch Dynamics) Model
}
\description{
Calculates the predicted catch, residuals, by fleet, population biomass, exploitation
rate, fishing mortality, of a fish stock using parameters estimated for a Catch
Dynamic Model.
}
\usage{
CatDynPred(x, method, partial = TRUE)
}
\arguments{
  \item{x}{
An object of class catdyn from function CatDynFit().
}
  \item{method}{
Character. The particular numerical method from which estimates have to be drawn.
}
  \item{partial}{
Logical, FALSE for transit models, TRUE (the default) for all other model versions.
}
}
\details{
This function is very similar to catdynexp() but instead of using arbitrary
parameter values given by the user, it takes the maximum likelihood estimates
produced by the call to the optimizer by function CatDynFit().
}
\value{
A list of length 2.
\item{Properties }{A list of length 3.
Units is a dataframe with the units of time step,
catch, body weight, and the numbers multiplier. Fleets is a dataframe with the
fleets names and the units of nominal effort for each fleet. Dates is a dataframe
with start and end dates of the fishing season in the ISO 8601 format.}
\item{Model }{A list of length 5.
Type is the perturbation type of model. Dates is the timing of perturbations,
Distr is the chosen likelihood model, Parameters is the
parameter values being explored, and Results, is a dataframe with the time step,
and for each fleet, observed effort, observed catch in numbers, predicted catch
in numbers, observed catch in weight, predicted catch in weight, deviance,
likelihood, deviance residuals, predicted population abundance, predicted population
biomass, observed exploitation rate, predicted exploitation rate, observed fishing
mortality, and predicted fishing mortality.}
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8), 1403-1415.

Roa-Ureta, R. H. et al. 2015. Fisheries Research 171 (Special Issue), 59-67.

Roa-Ureta, R. H. 2015. Fisheries Research 171 (Special Issue), 68-77.

Lin, Y-J. et al. 2017. Fisheries Research 195, 130-140.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
In resident stock models without emigration partial must be left at its default
value in order not to have any effect.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models  }