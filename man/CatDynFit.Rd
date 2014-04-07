\name{CatDynFit}
\alias{CatDynFit}
\title{
Fit Catch Dynamic Models by Maximum Likelihood
}
\description{
A wrapper and post-processing tool that checks that the data are passed with proper 
characteristics, calls optimx() (from package optimx) on any of 114 versions of the 
catch dynamics models (as internal functions), and then it post-processes optimx() 
results and join all results in a list of lists.
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
0 and 20. In two fleet cases 'p' is a two-components integer vector that quantifies
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
Character. Either 'normal' or 'lognormal' model for the observed random variable, 
the catch in numbers.
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
all model parameters. In multiannual applications and monthly time step (up to 20 
years of data) this might be time consuming but it should be carried out to
increase the chance that the optimizers will converge to resonable parameter space.

Initial parameter values must be passed log-transformed by the user. This is done to
be aware of scaling issues. CatDyn will backtransform the maximum likelihood estimates 
and its numerical Hessian matrix without user intervention using the delta method.

Generally, when 'p' is 5 (one fleet) or c(5,5) (two fleets) or lower, the model is 
applied to one annual season of data and the time step is day or week. Converserly, 
when 'p' is 6 (one fleet) or c(6,6) (two fleets) or higher the model is applied to 
multiannual series and the time step is the month, although it is conceivable that 
for a higly perturbed fishing system higher order 'p' values would be applied to single 
season cases.
}
\value{
A list of length 3.
\item{Data }{A list of length 2. 
'Properties' is a list of length 3. 'Units' is a dataframe with the units of time step, 
catch, body mass, and the numbers multiplier. 
'Fleets' is a dataframe with the fleets names and the units of nominal effort for 
each fleet. 
'Dates' is a dataframe with the start and end dates of the season in the ISO 8601 format. 
'Data' is a list of length equal to the number of fleets. Each component is a 
dataframe with the raw data, time step, observed effort, observed catch, observed 
mean body mass, observed catch in numbers, and the catch spike statistic.}
\item{Initial }{A dataframe with named initial values of all free parameters in the
model.}
\item{Model }{A list with length equal to the number of numerical methods. Each
component has the perturbation type model, the dates of events, the chosen distribution 
for the observation of catch, the integer code describing the succes or not of covergence 
returned by the method, the Karush Kuhn Tucker conditions, hopefully TRUE TRUE, the value 
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
optimization.
}
\examples{
#Falkland Islands squid fishery in 1990.
x <- lolgahi
#Note the second fleet data is fake
#Fit a 2-fleet 1P-1P model with lognormal observation error
#Create the data object
lgahi2F <- as.CatDynData(x=lolgahi,
                         step="day",
                         fleet.name=c("Fk","fake"),
                         coleff=c(2,5),
                         colcat=c(1,4),
                         colmbw=c(3,6),
                         unitseff=c("nboats","nfakes"),
                         unitscat="kg",
                         unitsmbw="kg",
                         nmult="bill",
                         season.dates=c(as.Date("1990-01-31"),
                                        as.Date("1990-05-30")))
#Interactive plot
#Examine the raw data
\dontrun{
plot(lgahi2F,mark=TRUE,offset=c(NA,NA,.75),hem="S")
}
#Set initial parameter values
M           <- 0.010 #1/Time step   
N0.ini      <- 6.5 #billions        
P1F1.ini    <- 1.2 #billions        
kF1.ini     <- 6.5e-05 #1/n of boats
alphaF1.ini <- 1.4 #adimensional    
betaF1.ini  <- 0.6 #adimensional    
P1F2.ini    <- 0.5 #billions        
kF2.ini     <- 5.5e-04 #1/n of boats
alphaF2.ini <- 0.4 #adimensional    
betaF2.ini  <- 0.8 #adimensional    
pars.ini  <- c(log(M),
               log(N0.ini),
               log(P1F1.ini),
               log(kF1.ini),
               log(alphaF1.ini),
               log(betaF1.ini),
               log(P1F2.ini),
               log(kF2.ini),
               log(alphaF2.ini),
               log(betaF2.ini))
#Set dates of events
P1F1  <- 70 #timing of perturbation of fleet 1
P1F2  <- 71 #timing of perturbation of fleet 2
dates <- c(head(lgahi2F$Data$Fk$time.step,1),
           P1F1,
           P1F2,
           tail(lgahi2F$Data$Fk$time.step,1))
#Create predictions from initial parameter values
lgahi2F.n.1P.ini <- catdynexp(x=lgahi2F,
                              p=c(1,1),
                              par=pars.ini,
                              dates=dates,
                              distr="lognormal") 
#Interactive plot
#Examine how good predictions are with the given initial parameter values
\dontrun{
plot(x=lgahi2F.n.1P.ini,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.4,
     Biom.ypos=0)
}
#Happy with initial parameter values, proceed to inference with two numerical
#methods
#Note, very low maximum number of iterations
(lgahi2F.ln.1P.fit <- CatDynFit(x=lgahi2F,
                                p=c(1,1),
                                par=pars.ini,
                                dates=dates,
                                distr="lognormal",
                                method=c("CG", "spg"),
                                itnmax=10))
#Examine results
lgahi2F.n.1P.pred.spg <- CatDynPred(lgahi2F.ln.1P.fit,"spg")
#Interactive plot
\dontrun{
plot(x=lgahi2F.n.1P.pred.spg,
     leg.pos="topright",
     Biom.tstep=7,
     Biom.xpos=0.18,
     Biom.ypos=0.1,
     AIC.xpos=0.18,
     AIC.ypos=0.2)
}
}
\keyword{ ~datasets }
\keyword{ ~manip }
\keyword{ ~optimize }
\keyword{ ~iplot }
\keyword{ ~distribution }
\keyword{ ~models }