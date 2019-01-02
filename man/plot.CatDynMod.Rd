\name{plot.CatDynMod}
\alias{plot.CatDynMod}
\title{
Examination of Results from Fitting Catch Dynamic Models
}
\description{
After model fit and prediction, examine model results on a graphical display.
}
\usage{
\method{plot}{CatDynMod}(x, leg.pos, Biom.tstep, Cat.tstep, Biom.xpos, Biom.ypos,
                         Cat.xpos, Cat.ypos,
                         diagnostics.panels = 2, top.lab = TRUE, ...)
}
\arguments{
  \item{x}{
An object of class CatDynMod.
}
  \item{leg.pos}{
The position of the legend in the top left panel. Passed to legend().
}
  \item{Biom.tstep}{
Integer. The number of time steps over which to average the population biomass 
counting from the end of the season backwards.
}
  \item{Cat.tstep}{
Integer. The number of time steps over which to add the catch counting from the
end of the season backwards.
}
  \item{Biom.xpos}{
Numeric. The position of the biomass estimate on the x-axis of the first panel 
of the plot in relative units.
}
  \item{Biom.ypos}{
Numeric. The position of the biomass estimate on the y-axis of the first panel 
of the plot in relative units.
}
  \item{Cat.xpos}{
Numeric. The position of the catch value on the x-axis of the first panel of
the plot in relative units.
}
  \item{Cat.ypos}{
Numeric. The position of the catch value on the y-axis of the first panel of
the plot in relative units.
}
  \item{diagnostics.panels}{
Integer. Whether to plot just one panel with the data and model prediction
(diagnostics.panels = 0), or add three more panels for diagnostics based on
deviance residual of the same size as the data-prediction panel
(diagnostics.panels = 1), or to make the deviance residual panels smaller
(diagnostics.panels = 2, the default).
}
  \item{top.lab}{
Logical. Whether to put a top label describing choices of model type, likelihood,
and numerical method. Defaults to TRUE.
}
  \item{\dots}{
Further arguments to pass to plot(), hist().
}
}
\details{
If the average population biomass over the whole season is to be posted then an
integer equal to the number of time steps in the season shall be entered for the
Biom.tstep argument.
}
\value{
A four panel plot of data, model predictions, and residual analysis.
}
\references{
Roa-Ureta, R. H. 2012. ICES Journal of Marine Science 69(8):1403-1415.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
The target symbols on the bottom of the tope left panel are the timings of any
perturbations set by the user.

In transit fisheries, entry target symbols are in red and exit target symbols 
are in blue.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~iplot }

