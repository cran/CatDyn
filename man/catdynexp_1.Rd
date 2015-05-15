\name{catdynexp}
\alias{catdynexp}
\title{
Exploratory Evaluation of Initial Values of CatDyn Models
}
\description{
Using a CatDynData object and initial parameter values, create predictions of
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
depletion. In one-fleet cases 'p' is a scalar integer that can take any value 
between -20 and 20. In two-fleet cases 'p' is a two-components integer vector that 
quantifies the number of perturbation events of each fleet. It can take values 
c(0,0), c(0,1), ..., c(0,5), c(1,1), ..., c(1,5), ..., c(4,5), c(5,5), c(6,6), 
..., c(20,20). In transit fisheries, where the stock only passes through the area 
where it is being fished and leave the fishing grounds at some time step within 
the season, 'p' should be negative and will take any integer value between -1 
and -20.
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
Character, either "poisson", "negbin", "normal","apnormal","lognormal","aplnormal",
or "gamma", or any pair of these seven (two-fleets systems), corresponding to the 
likelihood model.
}
}
\details{
The "negbin" value for the 'distr' parameter corresponds to the negative binomial
distribution for counts, as an alternative to the 'poisson' for cases where the
assumption of the mean equal to the variance is untenable.

The difference between "normal" and "apnormal", "lognormal" and "aplnormal" is 
that in the former the dispersion parameters is included in the likelihood function 
and it is a free parameter to be estimated along with the parameters of the 
generalized depletion model (and therefore an initial value for the dispersion must 
be provided) whereas in the latter the dispersion is eliminated by using the 
adjusted profile likelihood approximation. For the "negbin" and "gamma" distributions 
the dispersion parameter is always estimated along with the model parameters.
In two-fleets models any pair of the seven available likelihood models can be 
specified although since the units of the catch need to be the same for both fleets 
then it does not make sense to specify a combination of a distribution for counts 
and for a continuous variable.
}
\value{
A list of length 2.
\item{Properties }{A list of length 3. 
'Units' is a dataframe with the units of time step, 
catch, body weight, and the numbers multiplier. 'Fleets' is a dataframe with the 
fleets names and the units of nominal effort for each fleet. 'Dates' is a dataframe 
with start and end dates of the fishing season in the ISO 8601 format}
\item{Model }{A list of length 5.
'Type' is the perturbation type of model. 'Dates' is the timing of perturbations,
'Distr' is the chosen likelihood model, 'Parameters' is the
parameter values being explored, and 'Results', is a dataframe with the time step,
and for each fleet, the observed effort, observed catch, predicted catch, and the 
residuals, plus two more columns, one with the predicted population abundance, 
and a final one with the predicted population biomass.}
}
\author{
Ruben H. Roa-Ureta
}
\note{
The plot.CatDynData() function, acting on objects of class CatDynData, provides 
a plot of a statistic called the catch spike statistic, that can be useful to
determine the 'p' and 'dates' arguments.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models }

