\name{CatDynFit}
\alias{CatDynFit}
\title{
Fit CatDyn Models by Maximum Likelihood
}
\description{
A wrapper and post-processing tool that checks that the data are passed with proper 
characteristics, calls optimx() (from package optimx) on any of hundreds of possible
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
depletion. In one fleet cases 'p' is a scalar integer that can take any value between 
0 and 100. In two fleet cases 'p' is a two-components integer vector that quantifies
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
right after the time of entry.   For example, p=c(1,4,50,10,60,61) would specify
a two-perturbation model which starts at time step 1, has the first perturbation
at time step 4, with exit at time step 50, second perturbation at time 10, with
exit at time step 60, and season finishing at time step 61.
}
  \item{distr}{
Character, either 'normal', 'apnormal', 'lognormal', 'aplnormal', 'gamma', or 'poisson',
corresponding to the likelihood model.
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
all model parameters. In multiannual applications and monthly time step this might be
time consuming but it should be carried out to increase the chance that the optimizers
will converge to resonable parameter space.

Initial parameter values must be passed log-transformed by the user. CatDynFit() will
backtransform the maximum likelihood estimates and its numerical Hessian matrix
without user intervention using the delta method.

Generally, when 'p' is 5 (one fleet) or c(5,5) (two fleets) or lower, the model is 
applied to one annual season of data and the time step is day or week. Converserly, 
when 'p' is 6 (one fleet) or c(6,6) (two fleets) or higher the model is applied to 
multiannual series and the time step is the month, although it is conceivable that 
for a higly perturbed fishing system higher order 'p' values would be applied to single 
season cases.

The models set up for transit fisheries are single fleet only, so when 'p' is negative,
taking any value in the admissible range, its length must be 1.

The discrete Poisson distribution option is recommended for fisheries where the catch
is counted in number of fish instead of weight.
}
\value{
A list of length 3.
\item{Data }{A list of length 2. 
'Properties' is a list of length 3. 'Units' is a dataframe with the units of time step, 
catch, body weight, and the numbers multiplier.
'Fleets' is a dataframe with the fleets names and the units of nominal effort for 
each fleet. 
'Dates' is a dataframe with the start and end dates of the season in the ISO 8601 format. 
'Data' is a list of length equal to the number of fleets (1 or 2). Each component is a
dataframe with the raw data, time step, observed effort, observed catch, observed 
mean body weight, observed catch in numbers, and the catch spike statistic.}
\item{Initial }{A dataframe with named initial values of all free parameters in the
model.}
\item{Model }{A list with length equal to the number of numerical methods. Each
component has the perturbation type model, the dates of events, the chosen distribution 
for the observation of catch, the integer code describing the succes or not of covergence 
returned by the method, the Karush Kuhn Tucker conditions, hopefully TRUE and TRUE, the value
of the Akaike Information Criterion, not comparable between different distributions, the 
back-transformed (from log) maximum likelihood estimates, the numerical gradients at each 
maximum likelihood estimate, the standard errors of backtransformed (from log) maximum 
likelihood estimates, and the correlation matrix of the back-transformed (from log) maximum 
likelihood estimates.}
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
}
\author{
Ruben H. Roa-Ureta
}
\note{
Complex models may take several hours to converge on a PC. As an example, a
two fleet model with 18 perturbations each fleet, p=c(18,18), and a lognormal
observation model, totalling 44 parameters to estimate from 216 monthly time steps, 
coverged successfully in 16 hours on a Windows 7, 64 bit, 3 GHz processor, 8 GB RAM. 

Some effort has been made to avoid being kicked out of numerical optimization by
just one numerical method that fails, so that optimization continues with other 
methods, but there may remain some cases when the whole optimization process is 
aborted by failure in just one method. Try taking out some suspicious methods and 
optimize again. Experience shows that methods 'spg' and 'CG' are robust for this 
kind of model so both should be considered as the baseline case for numerical
optimization. When using the option of modeling transit fisheries with the Poisson
distribution it has been observed that methods 'bobyqa' and 'newuoa' outperform
'spg' and 'CG', so keep and open mind and take advantage of optimx by trying
several numerical methods.
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
                       season.dates=c(as.Date("1990-01-31"),
                                      as.Date("1990-05-30")))
plot(lgahi, mark=TRUE, offset=c(NA,NA,.75),hem="S")
#1 Fit a 1-fleet 1P model with normal observation error
M         <- 0.011   #1/Time step
N0.ini    <- 3.8     #billions
P1.ini    <- 1.3     #billions
k.ini     <- 5.0e-05 #1/n of boats
alpha.ini <- 1.7     #adimensional
beta.ini  <- 0.6     #adimensional
pars.ini  <- c(log(M),
               log(N0.ini),
               log(P1.ini),
               log(k.ini),
               log(alpha.ini),
               log(beta.ini))
#Dates
P1    <- 70 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
           P1,
           tail(lgahi$Data$Fk$time.step,1))
lgahi.n.1P.ini <- catdynexp(x=lgahi,
                            p=1,
                            par=pars.ini,
                            dates=dates,
                            distr="normal")
plot(x=lgahi.n.1P.ini,
     leg.pos="topright",
     Biom.tstep="mid",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit; note that the number of iterations, itnmax, should be much higher in
#real applications
(lgahi.n.1P.fit <- CatDynFit(x=lgahi,
                             p=1,
                             par=pars.ini,
                             dates=dates,
                             distr="normal",
                             method="spg",
                             itnmax=10))
#examine results
lgahi.n.1P.pred.spg <- CatDynPred(lgahi.n.1P.fit,"spg")
plot(x=lgahi.n.1P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#2 Fit a 1-fleet 1P model with normal observation error but the adjusted profile
#approximation
M         <- 0.012   #1/Time step
N0.ini    <- 4.3     #billions
P1.ini    <- 1.3     #billions
k.ini     <- 1.7e-04 #1/n of boats
alpha.ini <- 1.3     #adimensional
beta.ini  <- 0.6     #adimensional
pars.ini  <- c(log(M),
               log(N0.ini),
               log(P1.ini),
               log(k.ini),
               log(alpha.ini),
               log(beta.ini))
#Dates
P1    <- 70 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
          P1,
          tail(lgahi$Data$Fk$time.step,1))
lgahi.apn.1P.ini <- catdynexp(x=lgahi,
                                p=1,
                                par=pars.ini,
                                dates=dates,
                                distr="apnormal")
plot(x=lgahi.apn.1P.ini,
     leg.pos="topright",
     Biom.tstep="fin",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit; note that the number of iterations, itnmax, should be much higher in
#real applications
(lgahi.apn.1P.fit <- CatDynFit(x=lgahi,
                               p=1,
                               par=pars.ini,
                               dates=dates,
                               distr="apnormal",
                               method=c("spg", "CG"),
                               itnmax=10))
lgahi.apn.1P.pred.spg <- CatDynPred(lgahi.apn.1P.fit,"spg")
plot(x=lgahi.apn.1P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#3 Fit a 1-fleet 2P model with lognormal observation error with the adjusted
#profile approximation
M         <- 0.011   #1/Time step
N0.ini    <- 3.8     #billions
P1.ini    <- 1.3     #billions
P2.ini    <- 0.5     #billions
k.ini     <- 4.0e-05 #1/n of boats
alpha.ini <- 1.7     #adimensional
beta.ini  <- 0.6     #adimensional
pars.ini  <- c(log(M),
               log(N0.ini),
               log(P1.ini),
               log(P2.ini),
               log(k.ini),
               log(alpha.ini),
               log(beta.ini))
#Dates
P1    <- 70  #Selected by visual inspection of standard plot
P2    <- 135 #Selected by visual inspection of standard plot
dates <- c(head(lgahi$Data$Fk$time.step,1),
          P1,
          P2,
          tail(lgahi$Data$Fk$time.step,1))
lgahi.apln.2P.ini <- catdynexp(x=lgahi,
                               p=2,
                               par=pars.ini,
                               dates=dates,
                               distr="aplnormal")
plot(x=lgahi.apln.2P.ini,
     leg.pos="topright",
     Biom.tstep="ini",
     Biom.xpos=0.4,
     Biom.ypos=0)
#fit; note that the number of iterations, itnmax, should be much higher in
#real applications
(lgahi.apln.2P.fit <- CatDynFit(x=lgahi,
                                p=2,
                                par=pars.ini,
                                dates=dates,
                                distr="aplnormal",
                                method="spg",
                                itnmax=10))
#examine results
lgahi.apln.2P.pred.spg <- CatDynPred(lgahi.apln.2P.fit,"spg")
plot(x=lgahi.apln.2P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
#Summary table for model selection
(lgahi.sum <- CatDynSum(x=list(lgahi.n.1P.fit,
                               lgahi.apn.1P.fit,
                               lgahi.apln.2P.fit),
                         season=1990,
                         meths=c("spg","CG","spg")))
#Plot for correlations among parameter estimates
CatDynCor(x=list(lgahi.n.1P.fit,
                 lgahi.apn.1P.fit,
                 lgahi.apln.2P.fit),
          ttl=c("Normal 1P","APNormal 1P","APLognormal 2P"),
          meths=c("spg","CG","spg"),
          arr=c(3,1))
}
\keyword{ ~datasets }
\keyword{ ~manip }
\keyword{ ~optimize }
\keyword{ ~iplot }
\keyword{ ~distribution }
\keyword{ ~models }