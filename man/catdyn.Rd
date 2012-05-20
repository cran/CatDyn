\name{catdyn}
\alias{catdyn}
\title{
Fit Catch Dynamic Models by Maximum Likelihood
}
\description{
A wrapper that checks that the data are passed with proper dimensions, calls
optimx and any of 10 versions of the catch dynamics models, and then it
post-processes optimx or optim results and join all results in a list of lists.
}
\usage{
catdyn(p, par, dates, obseff, obscat, M.fixed, M, distr, num.wrapper, method, 
       control = list(), hessian, itnmax)
}
\arguments{
  \item{p}{
Integer, the number of positive perturbations, from 0 to 4.
}
  \item{par}{
Numeric vector of parameter estimates in the log scale.
}
  \item{dates}{
Integer vector with the time steps of start of season, perturbations (if any), and end of season.
}
  \item{obseff}{
A numeric vector with observed nominal effort by time step.
}
  \item{obscat}{
A numeric vector with observed catch in numbers by time step.
}
  \item{M.fixed}{
Logical, shall natural mortality be estimated (M.fixed=FALSE) or remain fixed at a given value?
}
  \item{M}{
A positive number or NULL; if M.fixed is TRUE, this number sets the value of the natural mortality rate.
}
  \item{distr}{
Character, either 'normal' or 'lognormal' model for the observed random variable, the catch.
}
  \item{num.wrapper}{
Character, either 'optimx' or 'optim', the wrapper for numerical optimization
}
  \item{method}{
Character vector. Any method accepted by optimx() can be used, but some may return warnings or errors.
}
  \item{control}{
A list of control arguments to be passed to optimx().
}
  \item{hessian}{
Logical. TRUE to estimate the hessian. Not implemented. Must be TRUE.
}
  \item{itnmax}{
Numeric, maximum number of iterations, to pass to optimx().
}
}
\details{
The function controls the estimation process. From a valid value of p, it will define the proper process model, either pure
depletion or any of 4 perturbation models, and will determine the distribution of the catch data, either normal or lognormal.
In the case of calling optimx with optimization methods based on gradients, it will use numerical gradients but it will also 
return analytical gradients evaluated at the maximum likelihood estimates. In the alternative case of calling optim the 
corresponding analytical gradient function has to be passed and the value of those gradients at the maximum likelihood estimates
is returned. The function will re-organize optimx's or optim's output by adding some items (AIC, back-transformed
MLEs, etc) and ignoring other items.
}
\value{
\item{model }{Type of model, 0, 1, 2, 3, or 4 perturbations}
\item{distr }{Distribition chosen, normal or lognormal}
\item{method }{Name of numerical optimization method}
\item{converg }{Convergence message}
\item{kkt }{The Karush-Kuhn-Tucker optimality conditions, if calling optimx}
\item{AIC }{The Akaike Information Criterion}
\item{bt.par }{Back-transformed maximum likelihood estimates of model parameters}
\item{num.grads }{Numerical gradients at the maximum likelihood estimates, if calling optimx}
\item{ana.grads }{Numerical gradients at the maximum likelihood estimates}
\item{bt.stdevs }{Estimated standard deviations of back-transformed maximum likelihood estimates of model parameters}
\item{Cor }{The estimated correlation matrix of maximum likelihood estimates of model parameters}
}
\author{
Ruben H. Roa-Ureta
}
\note{
The delta method from package msm is used to calculate back-transformed standard deviations.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Fit a 4-perturbation model with optimx - Using numerical gradients
#Initial parameter values
M                       <- 0.014 #1/Time step
N0.ini.B.1990.S1        <- 37.6 #billions
P1.ini.B.1990.S1        <- 0.25 #billions
P2.ini.B.1990.S1        <- 10.9 #billions
P3.ini.B.1990.S1        <- 0.01 #billions
P4.ini.B.1990.S1        <- 4.4 #billions
k.ini.B.1990.S1         <- 1e-05 #1/n of vessels
alpha.ini.B.1990.S1     <- 1.4 #adimensional
beta.ini.B.1990.S1      <- 0.9 #adimensional
pars.ini.B.1990.S1.4P   <- c(log(M),
                             log(N0.ini.B.1990.S1),
                             log(P1.ini.B.1990.S1),
                             log(P2.ini.B.1990.S1),
                             log(P3.ini.B.1990.S1),
                             log(P4.ini.B.1990.S1),
                             log(k.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
P1.B.1990.S1.4P         <- 56 #Selected by visual inspection of standard plot
P2.B.1990.S1.4P         <- 71
P3.B.1990.S1.4P         <- 86
P4.B.1990.S1.4P         <- 135
dates.B.1990.S1.4P      <- c(head(SeasonData.1990.S1.B$period,1),
                             P1.B.1990.S1.4P,
                             P2.B.1990.S1.4P,
                             P3.B.1990.S1.4P,
                             P4.B.1990.S1.4P,
                             tail(SeasonData.1990.S1.B$period,1))
#Try more numerical optimization methods, see help(optimx)
t.start                            <- Sys.time()
CatDynMod.1990.S1.4P.Normal.optimx <- catdyn(p=4,
                                             par=pars.ini.B.1990.S1.4P,
                                             dates=dates.B.1990.S1.4P,
                                             obscat=SeasonData.1990.S1.B$obscat,
                                             obseff=SeasonData.1990.S1.B$obseff2,
                                             M.fixed=FALSE,
                                             distr="normal",
                                             num.wrapper="optimx",                                             
                                             method=c("spg"),
                                             hessian=TRUE,
                                             itnmax=1000)
t.end                              <- Sys.time()
(t.process <- t.end-t.start)
#
CatDynMod.1990.S1.4P.Normal.optimx
#Fit with optim a 0-perturbation model - Using analytical gradients
#Try more numerical optimization methods, see help(optim)
M                       <- 0.01
N01.ini.B.1990.S1       <- 7.5
q.ini.B.1990.S1         <- 1e-5
alpha.ini.B.1990.S1     <- 1.3
beta.ini.B.1990.S1      <- 0.5
pars.ini.B.1990.S1.0P   <- c(log(M),
                             log(N01.ini.B.1990.S1),
                             log(q.ini.B.1990.S1),
                             log(alpha.ini.B.1990.S1),
                             log(beta.ini.B.1990.S1))
#Dates
dates.B.1990.S1.0P      <- c(head(SeasonData.1990.S1.B$period,1),
                             tail(SeasonData.1990.S1.B$period,1))
#Try more numerical optimization methods, see help(optimx)
t.start                            <- Sys.time()
CatDynMod.1990.S1.0P.LogNormal.optim  <- catdyn(p=0,
                                                par=pars.ini.B.1990.S1.0P,
                                                dates=dates.B.1990.S1.0P,
                                                obscat=SeasonData.1990.S1.B$obscat,
                                                obseff=SeasonData.1990.S1.B$obseff2,
                                                M.fixed=FALSE,
                                                distr="lognormal",
                                                num.wrapper="optim",
                                                method="BFGS",
                                                hessian=TRUE,
                                                control=list(maxit=1000),
                                                itnmax=NULL)
t.end                              <- Sys.time()
(t.process <- t.end-t.start)
#
CatDynMod.1990.S1.0P.LogNormal.optim


}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
