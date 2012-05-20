\name{plot.CatDynMod}
\alias{plot.CatDynMod}
\title{
Examining Model Predictions of Catch Dynamics Models on a Graphical Display
}
\description{
A four-panel plot showing the relation between catch model prediction and observed catch. The top-left panel shows catch data,
and model predictions along the season, plus the dates of perturbations, if any. It includes the AIC and the biomass at the start,
mid, or the end of the seasons. The top-right panel is a histogram of residuals. The bottom-left panel is a scatterplot of residuals.
The bottom-right panel is a Q-Q normal plot of residuals.
}
\usage{
  \method{plot}{CatDynMod}(x, tstep, mult, Biom, AIC, top.text, leg.pos, AIC.xpos,
               AIC.ypos, Biom.tstep, Biom.xpos, Biom.ypos, p.dates, \dots)
}
\arguments{
  \item{x}{
A data.frame, an object of class CatDynMod.
}
  \item{tstep}{
Character, the chosen time step (should be days or weeks).
}
  \item{mult}{
Character, the abundance multiplier.
}
  \item{Biom}{
Character, either the biomass at the start, at mid, or at the end of the season, or 'Collapse' (Abundance goes to zero).
}
  \item{AIC}{
The AIC of the model.
}
  \item{top.text}{
Character, a descriptive sentence.
}
  \item{leg.pos}{
Character, any of the choices for position of function legend.
}
  \item{AIC.xpos}{
Number, the position of the AIC value along the x-axis. 0 suppresses the text.
}
  \item{AIC.ypos}{
Number, the position of the AIC value along the y-axis. 0 suppresses the text.
}
  \item{Biom.tstep}{
Character, either "Ini" (or "ini"), "Mid" (or "mid"), or "Fin" (or "fin"), determining the timing of the biomass estimate.
}
  \item{Biom.xpos}{
Number, the position of the biomass value along the x-axis. 0 suppresses the text.
}
  \item{Biom.ypos}{
Number, the position of the biomass value along the y-axis. 0 suppresses the text.
}
  \item{p.dates}{
Integer or integer vector, the dates of perturbations. 0 for a pure depletion model.
}
  \item{\dots}{
Arguments to be passed to methods.
}
}
\details{
Argument x is a data.frame output from one of the functions CDMN0P, ..., CDMN4P.
}
\value{
A four-panel plot.
}
\author{
Ruben H. Roa-Ureta
}
\note{
This plot can be used to (1) choose good initial parameter values before optimization,
and to select time steps during the season when it is suspected that a perturbation might
have occurred (a large positive residual), and (2) for model checking after optimization.
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
#Parameters
M                       <- 0.014
N0.ini.B.1990.S1        <- 37.622
P1.ini.B.1990.S1        <- 0.25
P2.ini.B.1990.S1        <- 10.912
P3.ini.B.1990.S1        <- 0.010
P4.ini.B.1990.S1        <- 4.412
k.ini.B.1990.S1         <- 1.011e-5
alpha.ini.B.1990.S1     <- 1.433
beta.ini.B.1990.S1      <- 0.941
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
P1.B.1990.S1.4P         <- 56
P2.B.1990.S1.4P         <- 71
P3.B.1990.S1.4P         <- 86
P4.B.1990.S1.4P         <- 135
dates.B.1990.S1.4P      <- c(head(SeasonData.1990.S1.B$period,1),
                             P1.B.1990.S1.4P,
                             P2.B.1990.S1.4P,
                             P3.B.1990.S1.4P,
                             P4.B.1990.S1.4P,
                             tail(SeasonData.1990.S1.B$period,1))
sealen.B.1990.S1        <- dates.B.1990.S1.4P[2]-dates.B.1990.S1.4P[1]+1
##Catch Dynamics Matrix
B.1990.S1.4P.ini        <- CDMN4P(par=pars.ini.B.1990.S1.4P,
                                  dates=dates.B.1990.S1.4P,
                                  obscat=SeasonData.1990.S1.B$obscat,
                                  obseff=SeasonData.1990.S1.B$obseff2,
                                  obsmbm=SeasonData.1990.S1.B$obsmbm,
                                  M.fixed=FALSE,
                                  distr='normal')
#AIC
AIC.B.1990.S1.4P.ini    <- 2*length(pars.ini.B.1990.S1.4P) -
                           2*(-((sealen.B.1990.S1-2)/2)*
                            log(sum(B.1990.S1.4P.ini$resids^2)))
#Plot
plot(x=B.1990.S1.4P.ini,
     tstep='Day',
     mult='Billions',
     Biom=round(tail(B.1990.S1.4P.ini$npred,1)*1e9
                *mean(tail(SeasonData.1990.S1.B$obsmbm,7))*1e-3),
     AIC=round(AIC.B.1990.S1.4P.ini,1),
     top.text="Beauchene Season 1 1990 - 4P Model - Normal",
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(56,71,86,135))
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
