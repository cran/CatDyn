\name{catdyn}
\alias{catdyn}
\title{
Fit Catch Dynamic Models by Maximum Likelihood
}
\description{
A wrapper that check that the data are passed with proper dimensions, calls
optimx and any of 20 versions of the catch dynamic models, and then it
post-processes optimx results joining all results in a list of lists.
}
\usage{
catdyn(p, par, dates, obseff, obscat, M.fixed, M, distr, method,
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
Number, if M.fixed is FALSE; NULL, if it is TRUE
}
  \item{distr}{
Character, normal or lognormal model for the observed random variable, the catch.
}
  \item{method}{
Character vector. Any method accepted by optimx can be used.
}
  \item{control}{
A list of control arguments to be passed to optimx.
}
  \item{hessian}{
Logical. TRUE to estimate the hessian. Not implemented. Must be TRUE.
}
  \item{itnmax}{
Numeric, maximum number of iterations, to pass to optimx.
}
}
\details{
The function controls the estimation process. From a valid value of p, it will define the proper process model, either pure
depletion or any of 4 perturbation models, and will determine the distribution of the catch data, either normal or lognormal.
In the case of optimization methods based on gradients, it will use numerical gradients but it will also return analytical
gradients that are still under testing.

Note that some numerical optimization methods may return warnings or even errors, causing the function to stop. Try different
methods. Methods spg, CG, BFGS, Nelder-Mead, usually return results but these results nay differ substantially. It is
recommended to try different methods and compare results. It is also possible to use chains of methods. See help(optimx).

The function will re-organize optimx's output by adding some items (AIC, back-transformed
MLEs, etc) and ignoring other items.

The output is a list of results in a list of optimization methods.
}
\value{
\item{model }{Type of model, 0, 1, 2, 3, or 4 perturbations}
\item{distr }{Distribition chosen, normal or lognormal}
\item{method }{Name of numerical optimization method}
\item{M.fixed }{If M.fixed is TRUE this is the value of M}
\item{converg }{Convergence message}
\item{kkt }{The Karush-Kuhn-Tucker optimality conditions}
\item{AIC }{The Akaike Information Criterion}
\item{bt.par }{Back-transformed maximum likelihood estimates of model parameters}
\item{num.grads }{If gradient is 'numerical' this is the vector of numerical gradients at the maximum likelihood estimates, }
\item{ana.grads }{Analytical gradients at the maximum likelihood estimates}
\item{bt.stdevs }{Estimated standard deviations of back-transformed maximum likelihood estimates of model parameters}
\item{Cor }{The estimated correlation matrix of maximum likelihood estimates of model parameters}
}
\author{
Ruben H. Roa-Ureta
}
\note{
The deltamethod() function from package msm is used to calculate back-transformed standard deviations and correlation matrix.
This function may return warnings when the square root functions is applied on a NaN. Try increasing the number of iterations.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Example 1 - A Pure Depletion (0-Perturbation) Model
#Parameters - Initial par values found by using CDMN0P and plot.CatDynMod
M           <- 0.014 #1/Time step
N0          <- 7.5 #billions
k           <- 1.2e-5 #1/n of vessels
alpha       <- 1.32 #adimensional
beta        <- 0.45 #adimensional
pars.ini.0P <- c(log(M),log(N0),log(k),log(alpha),log(beta))
#Dates
dates.0P    <- c(head(SeasonData.1990.S1.B$period,1),
                 tail(SeasonData.1990.S1.B$period,1))
#Fit 1.1 - Normal Distribution with free M.
t.start             <- Sys.time()
CatDynMod.0P.Normal <- catdyn(p=0,
                              par=pars.ini.0P,
                              itnmax=500,
                              method=c("spg", "CG"),
                              hessian=TRUE,
                              dates=dates.0P,
                              obseff=SeasonData.1990.S1.B$obseff2,
                              obscat=SeasonData.1990.S1.B$obscat,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="normal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
CatDynMod.0P.Normal

#Fit 1.2 - LogNormal distribution with free M.
t.start                <- Sys.time()
CatDynMod.0P.LogNormal <- catdyn(p=0,
                                 par=pars.ini.0P,
                                 itnmax=500,
                                 method=c("spg", "CG"),
                                 hessian=TRUE,
                                 dates=dates.0P,
                                 obseff=SeasonData.1990.S1.B$obseff2,
                                 obscat=SeasonData.1990.S1.B$obscat,
                                 M.fixed=FALSE,
                                 M=NULL,
                                 distr="lognormal")
t.end                  <- Sys.time()
(t.process <- t.end-t.start)
CatDynMod.0P.LogNormal

#Fit 1.3 - LogNormal distribution with fixed M.
t.start                     <- Sys.time()
CatDynMod.0P.LogNormal.FixM <- catdyn(p=0,
                                      par=c(log(N0),
                                            log(k),
                                            log(alpha),
                                            log(beta)),
                                      itnmax=500,
                                      method=c("spg", "CG"),
                                      hessian=TRUE,
                                      dates=dates.0P,
                                      obseff=SeasonData.1990.S1.B$obseff2,
                                      obscat=SeasonData.1990.S1.B$obscat,
                                      M.fixed=TRUE,
                                      M=0.014,
                                      distr="lognormal")
t.end                       <- Sys.time()
(t.process <- t.end-t.start)
CatDynMod.0P.LogNormal.FixM

#Example 2 - A 1-Perturbation Model
#Parameters - Initial par values found by using CDMN0P and plot.CatDynMod
M           <- 0.009 #1/Time step
N0          <- 16.6 #billions
P1          <- 2.5 #billions
k           <- 1e-05 #1/n of vessels
alpha       <- 1.4 #adimensional
beta        <- 1.2 #adimensional
pars.ini.1P <- c(log(M),log(N0),log(P1),log(k),log(alpha),log(beta))
#Dates
P1.1P       <- 71 #Selected by visual inspection of standard plot
dates.1P    <- c(head(SeasonData.1990.S1.B$period,1),
                 P1.1P,
                 tail(SeasonData.1990.S1.B$period,1))
#Fit 2.1 - Normal distribution with fixed M.
t.start                     <- Sys.time()
CatDynMod.1P.Normal.FixM <- catdyn(p=1,
                                   par=c(log(N0),
                                         log(P1),
                                         log(k),
                                         log(alpha),
                                         log(beta)),
                                   method=c("spg", "CG"),
                                   hessian=TRUE,
                                   itnmax=500,
                                   dates=dates.1P,
                                   obscat=SeasonData.1990.S1.B$obscat,
                                   obseff=SeasonData.1990.S1.B$obseff2,
                                   M.fixed=TRUE,
                                   M=M,
                                   distr="normal")
t.end                    <- Sys.time()
(t.process <- t.end-t.start)
CatDynMod.1P.Normal.FixM

#Fit 2.2 - LogNormal distribution with free M.
t.start                        <- Sys.time()
CatDynMod.1P.LogNormal <- catdyn(p=1,
                                 par=pars.ini.1P,
                                 method=c("spg", "CG"),
                                 hessian=TRUE,
                                 itnmax=500,
                                 dates=dates.1P,
                                 obscat=SeasonData.1990.S1.B$obscat,
                                 obseff=SeasonData.1990.S1.B$obseff2,
                                 M.fixed=FALSE,
                                 M=NULL,
                                 distr="lognormal")
t.end                  <- Sys.time()
(t.process <- t.end-t.start)
CatDynMod.1P.LogNormal

}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
