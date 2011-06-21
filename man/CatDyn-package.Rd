\name{CatDyn-package}
\alias{CatDyn-package}
\alias{CatDyn}
\docType{package}
\title{
Fisheries Stock Assessment by Catch Dynamic Models
}
\description{
Using high-frequency catch and effort data (e.g. daily) CatDyn implements a type
of stock assessment models oriented to the operational fishing data. The estimated
parameters are in two groups, stock abundance and fishing operation. CatDyn includes
model versions that assume all the fish were available at the start of the fishing
season, or that up to 4 waves of recruitment or fleet area expansions occurred
during the season. CatDyn also includes two hypotheses for the distribution of the
observed response variable, either additive (normal) or multiplicative (lognormal).
The variance parameter in the chosen distribution is profiles out, and the profile
likelihood is adjusted. In CatDyn, the natural mortality rate parameter M can be
estimated with the other parameters or it can be fixed at a given value, making
inference on the other parameters conditional on the given fixed M value.
}
\details{
\tabular{ll}{
Package: \tab CatDyn\cr
Type: \tab Package\cr
Version: \tab 1.0-2\cr
Date: \tab 2011-05-27\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
Display exploratory plots to select a unit of effort, using
plot() with methods for class CatDynData.

Pick initial parameter values for one of the five versions of
the model and check them by using the exploratory plot produced by
plot() with methods for class CatDynMod.

Enter the selected initial parameter values transformed to the
log scale to the wrapper function catdyn() and fit the data to the model.

Inspect results using plot() on an object of class CatDynMod, using the
estimated parameter values.

The process equations in the catch Dynamic Models in this package are of the form:

C_t = k exp(-M/2) E_t^a N_t^b

N_t = N_0 exp(-Mt) + SUM_j P_j exp(-M(t-j)) - exp(-M/2) SUM_(i<t) C_(t-1) exp(-M(t-i-1))

for t=1, 2, 3, ...

where 'C' is catch in numbers, 't', 'i', 'j' are time step indicators, 'k' is a
scaling constant, 'E' is fishing effort, an observed predictor of catch, 'a' is a
parameter of effort synergy or saturability, 'N' is abundance, a latent predictor
of catch, 'b' is a parameter of hyperstability/hyperdepletion, and 'M' is natural
mortality rate in units of 1/time step. Latent abundance is expanded
into three processes. Initial abundance N_0 and its exponential decay due to
natural mortality M, a number (0, 1, 2, 3, or 4) of perturbations P_j that represent
fish migrations into the fishing grounds or expansions of the fishing
grounds by the fleet, resulting in point pulses of abundance, and a third term
that accounts for the fish that would have been dead due to natural causes but they
didn't have a chance because they were caught. This latter process is a discount
applied to the earlier catches in order to avoid an M-biased estimate of N_0.

}
\author{
Ruben H. Roa-Ureta

Maintainer: Ruben H. Roa-Ureta <ruben.roa.ureta@gmail.com>
}
\keyword{ package }
\examples{
#See examples in the function man pages, especially catdyn().
}
