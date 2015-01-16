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
0 and 100. In two-fleet cases 'p' is a two-components integer vector that quantifies
the number of perturbation events of each fleet. It can take values c(0,0), 
c(0,1), ..., c(0,5), c(1,1), ..., c(1,5), ..., c(4,5), c(5,5), c(6,6), ..., c(100,100).
In transit fisheries, where the stock only passes through the area where it is being
fished and leave the fishing grounds at some time step, 'p' should be negative and will
take any integer value between -20 and -1.
}
  \item{par}{
Numeric. Vector of initial parameter values in the log scale.
}
  \item{dates}{
Integer. Vector with the time steps of start of season, perturbations (if any), 
and end of season. In transit fisheries in addition to the timing of entry of
perturbations, the timing of exit for each perturbation shall also be provided,
right after the time of entry. For example, p=c(1,4,50,10,60,61) would specify
a two-perturbation model which starts at time step 1, has the first perturbation
at time step 4, with exit at time step 50, second perturbation at time 10, with
exit at time step 60, and season finishing at time step 61.
}
  \item{distr}{
Character, either 'normal','apnormal','lognormal','aplnormal','gamma','poisson',
corresponding to the likelihood model.
}
}
\details{
The plot.CatDynData() function, also acting on objects of class CatDynData, provides 
a plot of a statistic called the catch spike statistic, which can be useful to
determine the 'p' and 'dates' arguments.
The difference between 'normal' and 'apnormal' is that in the former the observation
variance is included in the likelihood function although it is not a free parameter
(it is computed after convergence from the estimated process parameters) whereas in
the latter the observation variance is eliminated by using the adjusted profile
likelihood approximation. The same applies to the difference between 'lognormal' and
'aplnormal',
}
\value{
A list of length 2.
\item{Properties }{A list of length 3. 
'Units' is a dataframe with the units of time step, 
catch, body weight, and the numbers multiplier. 'Fleets' is a dataframe with the fleets
names and the units of nominal effort for each fleet. 'Dates' is a dataframe with start 
and end dates of the fishing season in the ISO 8601 format}
\item{Model }{A list of length 5.
'Type' is the perturbation type of model. 'Dates' is the timing of perturbations,
'Distr' is the chosen likelihood model, 'Parameters' is the
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

