\name{plot.CatDynMod}
\alias{plot.CatDynMod}
\title{
Examination of Results from Fitting Catch Dynamic Models
}
\description{
After model fit and prediction, examine model results on a graphical display.
}
\usage{
\method{plot}{CatDynMod}(x, leg.pos, Biom.tstep, Biom.xpos, Biom.ypos, AIC.xpos, 
                         AIC.ypos, ...)
}
\arguments{
  \item{x}{
An object of class CatDynMod.
}
  \item{leg.pos}{
The position of the legend in the top left panel. Passed to legend().
}
  \item{Biom.tstep}{
Integer. The number of time steps over which to average the population biomass counting
from the end of the season or period backwards.
}
  \item{Biom.xpos}{
Numeric. The position of the biomass estimate on the x-axis of the first panel of the plot
in relative units.
}
  \item{Biom.ypos}{
Numeric. The position of the biomass estimate on the y-axis of the first panel of the plot
in relative units.
}
  \item{AIC.xpos}{
Numeric. The position of the Akaike Information Criterion on the x-axis of the first panel 
of the plot in relative units.
}
  \item{AIC.ypos}{
Numeric. The position of the Akaike Information Criterion on the y-axis of the first panel 
of the plot in relative units.
}
  \item{\dots}{
Further arguments to pass to plot(), hist().
}
}
\details{
If the average population biomass over the whole season is to be posted then an
integer equal to the number of time steps in the season shall be entered for the
'Biom.tstep' argument.
}
\value{
A four panel plot of data, model predictions, and residual analysis.
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
}
\author{
Ruben H. Roa-Ureta
}
\note{
The target symbols on the bottom of the first panel are the timings of any
perturbations set by the user.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~iplot }

