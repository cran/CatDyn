\name{CatDynFit}
\alias{CatDynFit}
\title{
Fit CatDyn Models by Maximum Likelihood
}
\description{
A wrapper and post-processing tool that checks that the data are passed with proper 
characteristics, calls optimx() (from package optimx) on any of dozens of possible
versions of the generalized depletion models (as internal functions), and then it
post-processes optimx() results and join all results in a list of lists.
}
\usage{
CatDynFit(x, p, par, dates, distr, method, control = list(), 
          hessian = TRUE, itnmax)
}
\arguments{
  \item{x}{
A data object of class CatDynData. See as.CatDynData().
}
  \item{p}{
Integer. The process model type, which quantifies the number of perturbations to
depletion. In one-fleet cases 'p' is a scalar integer that can take any value between
-20 and 20. In two-fleet cases 'p' is a two-components integer vector that quantifies
the number of perturbation events of each fleet. It can take values c(0,0),
c(0,1), ..., c(0,5), c(1,1), ..., c(1,5), ..., c(4,5), c(5,5), c(6,6), ..., c(20,20).
In transit fisheries, where the stock only passes through the area where it is being
fished and leave the fishing grounds at some time step within the season, 'p' should 
be negative and will take any integer value between -1 and -20.
}
  \item{par}{
Numeric. Vector of initial parameter values in the log scale.
}
  \item{dates}{
Integer. Vector with the time steps of start of season, perturbations (if any), 
and end of season. In transit fisheries in addition to the timing of entry of
perturbations, the timing of exit for each perturbation shall also be provided,
right after the time of entry.   For example, p=c(1,4,50,10,60,61) would specify
a two-perturbation model which starts at time step 1, has the first perturbation
at time step 4, with exit at time step 50, second perturbation at time 10, with
exit at time step 60, and season finishing at time step 61.
}
  \item{distr}{
Character, either "poisson", "negbin", "normal","apnormal","lognormal","aplnormal",
or "gamma" (one fleet), or any pair of these seven (two fleets) corresponding to 
the likelihood model.
}
  \item{method}{
Character. Any method accepted by optimx() can be used, but some may return 
warnings or errors.
}
  \item{control}{
A list of control arguments to be passed to optimx().
}
  \item{hessian}{
Logical. Defaults to TRUE. If set to FALSE all numerical methods tried will fail.
}
  \item{itnmax}{
Numeric. Maximum number of iterations, to pass to optimx().
}
}
\details{
Much care should be taken in selecting good initial values to pass in the par 
argument. To accomplish this CatDyn includes the CatDynExp class, and the catdynexp()
and the plot.CatDynExp() functions to graphically fine tune the initial values for
all model parameters. In multi-annual applications and monthly time step this might 
be time consuming but it should be carried out to increase the chance that the 
optimizers will converge to reasonable parameter space.

Initial parameter values must be passed log-transformed by the user. CatDynFit() will
backtransform the maximum likelihood estimates and its numerical Hessian matrix
without user intervention using the delta method.

Generally, when 'p' is 5 or lower (one fleet) or c(5,5) (two fleets) or lower, 
the model is applied to one annual season of data and the time step is "day" or 
"week". Conversely, when 'p' is 6 (one fleet) or c(6,6) (two fleets) or higher 
the model is applied to multiannual series and the time step is the month, although 
it is conceivable that for a higly perturbed fishing system higher 'p' values would 
be applied to single season cases.

The models set up for transit fisheries are single fleet only, so when 'p' is 
negative, taking any value in the admissible range, its length must be 1.

The discrete Poisson distribution option is recommended for fisheries where the 
catch is counted in number of fish instead of weight.

The "negbin" value for the 'distr' parameter corresponds to the negative binomial
distribution for counts, as an alternative to "poisson" for cases where the
assumption of the mean equal to the variance is untenable.

The difference between "normal" and "apnormal", "lognormal" and "aplnormal" is that
in the former the dispersion parameters is included in the likelihood function and 
it is a free parameter to be estimated along with the parameters of the generalized 
depletion model (and therefore an initial value for the dispersion has to be 
provided) whereas in the latter the dispersion is eliminated by using the adjusted 
profile likelihood approximation. For the "negbin" and "gamma" distributions the 
dispersion parameter is always estimated along with the model parameters.

In two-fleets models any pair of the seven available likelihood models can be 
specified although since the units of the catch need to be the same for both fleets 
then it does not make sense to specify a combination of a distribution for the 
catch in counts ("poisson" and "negbin") and a continuous distribution for the 
catch in biomass.
}
\value{
A list of length 3.
\item{Data }{A list of length 2. 
'Properties' is a list of length 3. 'Units' is a dataframe with the units of time 
step, catch, body weight, and the numbers multiplier.
'Fleets' is a dataframe with the fleets names and the units of nominal effort for 
each fleet. 
'Dates' is a dataframe with the start and end dates of the season in the ISO 8601 
format. 'Data' is a list of length equal to the number of fleets (1 or 2). Each 
component is a dataframe with the raw data, time step, observed effort, observed 
catch, observed mean body weight, observed catch in numbers, and the catch spike 
statistic.}
\item{Initial }{A dataframe with named initial values of all free parameters in 
the model.}
\item{Model }{A list with length equal to the number of numerical methods. Each
component has the perturbation type model, the dates of events, the chosen 
distribution for the observation of catch, the integer code describing the success 
or not of covergence returned by the method, the Karush Kuhn Tucker conditions, 
hopefully TRUE and TRUE, the value of the Akaike Information Criterion, not 
comparable between different distributions, the back-transformed (from log) maximum 
likelihood estimates, the numerical gradients at each maximum likelihood estimate, 
the standard errors of backtransformed (from log) maximum 
likelihood estimates, and the correlation matrix of the back-transformed (from log) 
maximum likelihood estimates.}
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
Roa-Ureta, R. H. 2015. Fisheries Research (In Press), doi:10.1016/j.fishres.2014.08.014
Roa-Ureta, R. H. et al. 2015. Fisheries Research (In Press), doi:10.1016/j.fishres.2014.12.006
}
\author{
Ruben H. Roa-Ureta
}
\note{
Complex models may take several hours to converge on a PC. As an example, a
two fleet model with 18 perturbations each fleet, p=c(18,18), and the 'aplnormal'
likelihood model, totalling 44 parameters to estimate from 216 monthly time 
steps, coverged successfully in 16 hours on a Windows 7, 64 bit, 3 GHz processor, 
8 GB RAM. 

Some effort has been made to avoid being kicked out of numerical optimization by
just one numerical method that fails, so that optimization continues with other 
methods, but there may remain some cases when the whole optimization process is 
aborted by failure in just one method. Try taking out some suspicious methods and 
optimize again. Experience shows that methods 'spg' and 'CG' are robust for this 
kind of model so both should be considered as the baseline case for numerical
optimization. When using the option of modeling transit fisheries with the Poisson
distribution it has been observed that methods 'bobyqa' and 'newuoa' also perform
well, so keep an open mind and take advantage of optimx by trying several numerical
methods.
}
\examples{
#Falkland Islands one-fleet squid fishery in 1990.
#Create the data object
lgahi <- as.CatDynData(x=lolgahi,
                       step="day",
                       fleet.name="Fk",
                       coleff=2,
                       colcat=1,
                       colmbw=3,
                       unitseff="nboats",
                       unitscat="kg",
                       unitsmbw="kg",
                       nmult="bill",
                       season.dates=c(as.Date("1990-01-31"),as.Date("1990-05-30")))
plot(lgahi,mark=TRUE,offset=c(NA,NA,.75),hem="S")
#
#1) Fit a 1-fleet 1P model with lognormal observation error and the adjusted
#profile approximation to the likelihood to eliminate the dispersion parameter
M         <- 0.011 #1/Time step
N0.ini    <- 3.8 #billions
P1.ini    <- 1.3 #billions
k.ini     <- 5.0e-05 #1/n of boats
alpha.ini <- 1.7 #adimensional
beta.ini  <- 0.6 #adimensional
pars.ini  <- log(c(M,
                   N0.ini,
                   P1.ini,
                   k.ini,
                   alpha.ini,
                   beta.ini))
#Dates
P1    <- 70 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
           P1,
           tail(lgahi$Data$Fk$time.step,1))
lgahi.apln.1P.ini <- catdynexp(x=lgahi,
                               p=1,
                               par=pars.ini,
                               dates=dates,
                               distr="aplnormal")
plot(x=lgahi.apln.1P.ini,
     leg.pos="topright",
     Biom.tstep="mid",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit
(lgahi.apln.1P.fit <- CatDynFit(x=lgahi,
                                p=1,
                                par=pars.ini,
                                dates=dates,
                                distr="aplnormal",
                                method="spg",
                                itnmax=10))
#examine results
lgahi.apln.1P.pred.spg <- CatDynPred(lgahi.apln.1P.fit,"spg")
plot(x=lgahi.apln.1P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#
#2) Fit a 1-fleet 2P model with lognormal observation error and full exact
#likelihood including the dispersion parameter
M         <- 0.011 #1/Time step
N0.ini    <- 3.8 #billions
P1.ini    <- 1.3 #billions
P2.ini    <- 0.5 #billions
k.ini     <- 4.0e-05 #1/n of boats
alpha.ini <- 1.7 #adimensional
beta.ini  <- 0.6 #adimensional
#Note how to get reasonable initial value for dispersion parameter
psi.ini   <- 0.33*sd(log(lgahi$Data$Fk$obscat.bill))^2
pars.ini  <- log(c(M,
                   N0.ini,
                   P1.ini,
                   P2.ini,
                   k.ini,
                   alpha.ini,
                   beta.ini,
                   psi.ini))
#Dates
P1    <- 70  #Selected by visual inspection of standard plot
P2    <- 135 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
           P1,
           P2,
           tail(lgahi$Data$Fk$time.step,1))
lgahi.ln.2P.ini <- catdynexp(x=lgahi,
                             p=2,
                             par=pars.ini,
                             dates=dates,
                             distr="lognormal")
plot(x=lgahi.ln.2P.ini,
     leg.pos="topright",
     Biom.tstep="ini",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit lognormal
(lgahi.ln.2P.fit <- CatDynFit(x=lgahi,
                              p=2,
                              par=pars.ini,
                              dates=dates,
                              distr="lognormal",
                              method="spg",
                              itnmax=10))
#examine results
lgahi.ln.2P.pred.spg <- CatDynPred(lgahi.ln.2P.fit,"spg")
plot(x=lgahi.ln.2P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#
#3) Fit a 1-fleet 2P model with gamma observation error and full, exact
#likelihod including the dispersion parameter
M         <- 0.011 #1/Time step
N0.ini    <- 3.8 #billions
P1.ini    <- 1.3 #billions
P2.ini    <- 0.5 #billions
k.ini     <- 4.0e-05 #1/n of boats
alpha.ini <- 1.7 #adimensional
psi.ini   <- 0.08 #From psi max. lik. estimate in 2P lognormal model
pars.ini  <- log(c(M,
                   N0.ini,
                   P1.ini,
                   P2.ini,
                   k.ini,
                   alpha.ini,
                   beta.ini,
                   psi.ini))
#Dates
P1    <- 70  #Selected by visual inspection of standard plot
P2    <- 135 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
           P1,
           P2,
           tail(lgahi$Data$Fk$time.step,1))
lgahi.g.2P.ini <- catdynexp(x=lgahi,
                            p=2,
                            par=pars.ini,
                            dates=dates,
                            distr="gamma")
plot(x=lgahi.g.2P.ini,
     leg.pos="topright",
     Biom.tstep="ini",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit gamma
(lgahi.g.2P.fit <- CatDynFit(x=lgahi,
                             p=2,
                             par=pars.ini,
                             dates=dates,
                             distr="gamma",
                             method="spg",
                             itnmax=10))
#examine results
lgahi.g.2P.pred.spg <- CatDynPred(lgahi.g.2P.fit,"spg")
plot(x=lgahi.g.2P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#
#Summary table for model selection
(lgahi.sum <- CatDynSum(x=list(lgahi.apln.1P.fit,
                               lgahi.ln.2P.fit,
                               lgahi.g.2P.fit),
                         season=1990,
                         meths=c("spg","spg","spg")))
#Plot for correlations among parameter estimates
#
CatDynCor(x=list(lgahi.apln.1P.fit,
                 lgahi.ln.2P.fit,
                 lgahi.g.2P.fit),
          ttl=c("Adjusted Profile Lognormal 1P","Lognormal 2P","Gamma 2P"),
          meths=c("spg","spg","spg"),
          arr=c(2,2))
}
\keyword{ ~datasets }
\keyword{ ~manip }
\keyword{ ~optimize }
\keyword{ ~iplot }
\keyword{ ~distribution }
\keyword{ ~models }