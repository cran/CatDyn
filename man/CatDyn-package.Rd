\name{CatDyn-package}
\alias{CatDyn-package}
\docType{package}
\title{
Fisheries Stock Assessment by Generalized Depletion Models
}
\description{
Using high-frequency (daily, weekly) or medium frequency (monthly) catch and effort 
data CatDyn implements a type of stock assessment model oriented to the operational 
fishing data. The estimated parameters are in two groups, stock abundance and fishing
operation. CatDyn includes 237 versions of the models depending on the number
of fleets, the number of perturbations to depletion, whether the stock is resident or
in transit, and six likelihood models.
}
\details{
\tabular{ll}{
Package: \tab CatDyn\cr
Type: \tab Package\cr
Version: \tab 1.0-5\cr
Date: \tab 2014-04-06\cr
License: \tab GPL (>= 2)\cr
}
Create a data object using raw data and the as.CatDynData() function.
Examine the data for regularities and perturbations using the generic plot()
function on an object of class CatDynData.
Examine the goodness of initial parameter values before statistical inference
by using the catdynexp() exploratory prediction function and the plot
generic function on an object of class CatDynExp.
Fit the model to the data by using the wrapper function CatDynFit(), which in turn
will call the optimx() optimizer wrapper of package optimx, with several 
numerical methods available to be used.
Examine the quality of the fit with the plot() function on an object of class
CatDynMod created by the CatDynPred() function.
Compare different models fit to the data with CatDynSum() and CatDynCor, based
on information theoretic, statistical, and numerical criteria.

The process equations in the Catch Dynamics Models in this package are of the form
\deqn{C_t = k e^{-M/2} E^{a}_{t} N^{b}_{t}}
\deqn{N_{t}=N_0 e^{-Mt}-e^{M/2}\sum_{i<t} C_{t-1} e^{-M(t-i-1)}+\sum_{j}P_{j} e^{-M(t-j)}}
where \emph{C} is catch in numbers, \emph{t}, \emph{i} are time step indicators, \emph{j} 
is perturbation index (\emph{j}=1,2,...,100), \emph{k} is a scaling constant, \emph{E} is nominal
fishing effort, an observed predictor of catch, \emph{a} is a parameter of effort synergy 
or saturability, \emph{N} is abundance, a latent predictor of catch, \emph{b} is a parameter 
of hyperstability or hyperdepletion, and \emph{M} is natural mortality rate per time step.
The second summand of the expanded latent predictor is a discount
applied to the earlier catches in order to avoid an \emph{M}-biased estimate of
initial abundance.
Perturbations to depletion represent fish migrations into the fishing grounds or expansions of the fishing
grounds by the fleet(s) resulting in point pulses of abundance.
In transit models (limited to up to 20 perturbations and one fleet) there are also emmigration events happening at
specific time steps for each perturbation.
In 2 fleet cases the fleets contribute complementary information about stock abundance, and
thus operate additively; any interaction between the fleets is latent and affects the estimated values
of fleet dependent parameters, such as \emph{k}, \emph{a}, and \emph{b}.

The observation model can take any of the following forms: an additive random normal term
added to the catch predicted by the process (normal and adjusted profile normal), a multiplicative
exponential term acting on the process-predicted catch such as the logarithm
of this multiplier distributes normally (lognormal and adjusted profile lognornmal), Gamma (shape and
scale parameterization), or a basic count process yielding a Poisson distribution.

}
\author{
Ruben H. Roa-Ureta <ruben.roa.ureta@mail.com>
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
}
\keyword{ package }
\keyword{ optimize}
\keyword{ models}
\keyword{ nonlinear}
\keyword{ datasets}
\keyword{ manip}
\examples{
#See examples for CatDynFit()
}
