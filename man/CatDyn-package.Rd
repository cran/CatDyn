\name{CatDyn-package}
\alias{CatDyn-package}
\docType{package}
\title{
Fisheries Stock Assessment by Generalized Depletion (Catch Dynamics) Models
}
\description{
Using high-frequency (daily, weekly) or medium frequency (monthly) catch and effort 
data CatDyn implements a type of stock assessment model oriented to the operational 
fishing data. The estimated parameters are in two groups, stock abundance (initial
abundance, episodic pulses of abundance, and natural mortality) and fishing operation
(hyperstability/hyperdepletion, saturability/synergy, catchability). CatDyn includes
100 versions of the models depending on the number of fleets (1 or 2), the number of
perturbations to depletion (1 to 25), whether the stock is resident without emigration,
resident with emigration, or in transit, and 72 likelihood models for the data depending
on the kind of random mechanism(s) assumed to have generated the data. The package has
graphical functions for exploratory analysis (to select reasonable initial values of
parameter to me estimated), it calls the optimx function from package optimx to use
many options for the numerical method to minimize the negative loglikelihood, and
several functions for post-analysis of results to select a best model for the data.
Considering the combinations that can be made with model versions, likelihood options,
and numerical methods, any data set can be fit to hundreds of model varieties making
model selection a substantial part of the stock assessment work.

Once the best model has been selected, the package also has functions to produce the
input information to fit a population dynamics model in a hierarchical inference
framework. See Roa-Ureta et al. (2015) below.
}
\details{
\tabular{ll}{
Package: \tab CatDyn\cr
Type: \tab Package\cr
Version: \tab 1.1-1\cr
Date: \tab 2018-12-20\cr
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
Compare different models fit to the data with CatDynSum(), CatDynCor(), and
CatDynPar(), based on information theoretic, statistical, and numerical criteria.
Create the input file for further population dynamics modeling with CatDynBSD().

The process equations in the Catch Dynamics Models in this package are of the form
\deqn{C_t = k e^{-M/2} E^{a}_{t} N^{b}_{t}}
whith versions for resident stock without emigration
\deqn{N_{t}=N_0 e^{-Mt}-e^{M/2}\sum_{i<t} C_{t-1} e^{-M(t-i-1)}+\sum_{j}R_{j} e^{-M(t-j)}}
resident stock with emigration
\deqn{N_{t}=N_0 e^{-Mt}-e^{M/2}\sum_{i<t} C_{t-1} e^{-M(t-i-1)}+\sum_{j}R_{j} e^{-M(t-v_j)}-\sum_{l}S_{l} e^{-M(t-w_l)}}
and stock in transit through the fishing grounds
\deqn{N_{t}=N_0 e^{-Mt}-e^{M/2}\sum_{i<t} C_{t-1} e^{-M(t-i-1)}+\sum_{j}I_j R_{j} e^{-M(t-v_j)}-\sum_{l}J_j R_{j} e^{-M(t-w_j)}}
where \emph{C} is the observed catch in numbers, a random variable, \emph{t},
\emph{i}, \emph{l}, are time step indicators, \emph{j}
is perturbation index (\emph{j}=1,2,...,25), \emph{k} is a scaling constant, \emph{E} is
nominal fishing effort, a predictor of catch observed exactly, \emph{a} is a parameter of effort
synergy or saturability, \emph{N} is abundance, a latent predictor of catch, \emph{b} is
a parameter of hyperstability or hyperdepletion, \emph{M} is natural mortality rate
per time step, \emph{I} and \emph{J} are indicator variables, taking values of 0 before
\emph{v} and 1 after, before \emph{w} and 1 after, respectively. The second summand of the expanded
latent predictor is a discount applied to the earlier catches in order to avoid an
\emph{M}-biased estimate of initial abundance.

Positive perturbations to depletion represent fish migrations into the fishing grounds or expansions
of the fishing grounds by the fleet(s) resulting in point pulses of abundance. Negative perturbations
represente emigration of parts of the stock (resident stock with emigration above) or the totality
of a recruitment wave (transit stock above).

In 2 fleet cases the fleets contribute complementary information about stock abundance, and
thus operate additively; any interaction between the fleets is latent and affects the 
estimated values of fleet dependent parameters, such as \emph{k}, \emph{a}, and \emph{b}.

The catch observation model can take any of the following forms: a Poisson counts process or a
negative binomial counts process for catch recorded in numbers, an additive random normal 
term added to the continuous catch (in weight) predicted by the process (normal and adjusted 
profile normal), a multiplicative exponential term acting on the process-predicted catch such 
as the logarithm of this multiplier distributes normally (lognormal, adjusted profile
lognornmal, and robust lognormal), Gamma (shape and scale parameterization), and Gumbel.
}
\author{
Ruben H. Roa-Ureta <ruben.roa.ureta@mail.com> (ORCID ID 0000-0002-9620-5224)
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8), 1403-1415.

Roa-Ureta, R. H. et al. 2015. Fisheries Research 171 (Special Issue), 59-67.

Roa-Ureta, R. H. 2015. Fisheries Research 171 (Special Issue), 68-77.

Lin, Y-J. et al. 2017. Fisheries Research 195, 130-140.
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
