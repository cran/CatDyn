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
catdynexp(x, p, par, dates, distr, partial=TRUE)
}
\arguments{
  \item{x}{
An object of class CatDynData.
}
  \item{p}{
Integer. The process model type, which quantifies the number of perturbations to
depletion. In one-fleet cases p is a scalar integer that can take any value
between -25 and 25. In two-fleet cases p is a two-components integer vector that
quantifies the number of perturbation events of each fleet. It can take values
c(0,0), c(0,1), ..., c(0,5), c(1,1), ..., c(1,5), ..., c(4,5), c(5,5), c(6,6),
..., c(25,25). In fisheries with emigration, where in addition to perturbations due
to positive pulses of abundance, there are perturbation due to negative pulses of
exodus, p should be negative and will take any integer value between -1 and -25,
this number fixing both the number of input and exit pulses.
}
  \item{par}{
Numeric. Vector of initial parameter values in the log scale.
}
  \item{dates}{
Integer. Vector with the time steps of start of season, perturbations (if any), 
and end of season. In fisheries with emigration, in addition to the timing of entry
of perturbations, the timing of exit for each perturbation shall also be provided,
right after the time of entry. For example, p=c(1,4,50,10,60,61) would specify
a two-perturbations model which starts at time step 1, has the first input
perturbation at time step 4, first exit perturbation at time step 50, second
input perturbation at time 10, second exit perturbation at time step 60, and
season finishing at time step 61.
}
  \item{distr}{
Character, either "poisson", "negbin", "normal", "apnormal", "lognormal",
"aplnormal", "gamma", "roblognormal", "gumbel", or any pair of these seven
(2-fleets systems), corresponding to the likelihood model.
}
  \item{partial}{
Logical, if the model incluyes emigration (p between -1 and -25) partial = TRUE
(the default) means that an unknown part of the abundance emigrate at some time
step within the season, causing a negative pulse of abundance. Alternatively,
partial = FALSE is used for transit stock models where the whole pulse of input
yet surviving exits the fishing grounds.
}
}
\details{
The "negbin" value for the distr parameter corresponds to the negative binomial
distribution for counts, as an alternative to the poisson for cases where the
assumption of the mean equal to the variance is untenable.

The difference between "normal" and "apnormal", "lognormal" and "aplnormal" is 
that in the former the dispersion parameters is included in the likelihood function 
and it is a free parameter to be estimated along with the parameters of the 
generalized depletion model (and therefore an initial value for the dispersion must 
be provided) whereas in the latter the dispersion is eliminated by using the 
adjusted profile likelihood approximation. Setting distr="roblognormal" will
use a robustified version of the lognormal which includes the dispersion parameter.
For the "negbin", "gamma" and "gumbel" distributions the dispersion parameter is
always estimated along with the model parameters. In two-fleets models any pair of
the nine available likelihood models can be specified.

In models with emigration, the logical parameter partial should be TRUE (the default)
when the magnitude of the exit pulse is an unknown parameter that needs to be estimated.
This could be the case of females exiting the fishing grounds to spawn, for instance.
On the other hand partial should be FALSE in transit stock fisheries, and in this case
the magnitude of the exit pulse need not be estimated, since all survivors leave.
Thus setting partial to TRUE and p negative estimates more parameters than
non-emigration and transit stock models.

At this point there is no allowance for models with emigration that include a different
number of entry and exit pulses. Transit stock fisheries (partial=FALSE) must have the
same number of entry and exit pulses but other emigration models (partial=TRUE) need not
be restricted. Next versions of the function will include options for different
number of entry and exit pulses.
}
\value{
A list of length 2.
\item{Properties }{A list of length 3. 
Units is a dataframe with the units of time step,
catch, body weight, and the numbers multiplier. Fleets is a dataframe with the
fleets names and the units of nominal effort for each fleet. Dates is a dataframe
with start and end dates of the fishing season in the ISO 8601 format}
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
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
The plot.CatDynData() function, acting on objects of class CatDynData, provides 
a plot of a statistic called the catch spike statistic, that can be useful to
determine the p and dates arguments.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~models }

