\name{catdynexp}
\alias{catdynexp}
\title{
Exploratory Evaluation of Initial Values of CatDyn Models
}
\description{
Using a CatDynData object and initial parameter values create predictions of
model results to examine the goodness of initial parameter values before
passing them to the numerical optimizer.
}
\usage{
catdynexp(x, p, par, dates, distr)
}
\arguments{
  \item{x}{
An object of class CatDynData.
}
  \item{p}{
Integer. The process model type, which quantifies the number of perturbations to
depletion. In one-fleet cases 'p' is a scalar integer that can take any value between 
0 and 20. In two-fleet cases 'p' is a two-components integer vector that quantifies
the number of perturbation events of each fleet. It can take values c(0,0), 
c(0,1), ..., c(0,5), c(1,1), ..., c(1,5), ..., c(4,5), c(5,5), c(6,6), ..., c(20,20).
}
  \item{par}{
Numeric. Vector of initial parameter values in the log scale.
}
  \item{dates}{
Integer. Vector with the time steps of start of season, perturbations (if any), 
and end of season.
}
  \item{distr}{
Character, either 'normal' or 'lognormal' model for the observed random variable, 
the catch.
}
}
\details{
The plot.CatDynData() function, also acting on objects of class CatDynData, provides 
a plot of a statistic called the catch spike statistic, which can be useful to
determine the 'p' and 'dates' arguments.
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
'Distr' is the chosen distribution for the observed catch data, 'Parameters' is the 
parameter values being explored, and 'Results', is a dataframe with the time step, 
and for each fleet, the observed effort, observed catch, predicted catch, and the 
residuals, plus two more columns, one with the predicted population abundance, 
and a final one with the predicted population biomass. }
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models }

