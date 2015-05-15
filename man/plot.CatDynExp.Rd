\name{plot.CatDynExp}
\alias{plot.CatDynExp}
\title{
Exploratory Modeling of Fisheries with Catch Dynamic Models
}
\description{
Evaluate and refine the goodness of initial parameter values before fitting 
catch dynamic models to data.
}
\usage{
\method{plot}{CatDynExp}(x, leg.pos, Biom.tstep, Biom.xpos, Biom.ypos, ...)
}
\arguments{
  \item{x}{
An object of class CatDynExp
}
  \item{leg.pos}{
The position of the legend in the first panel. Passed to legend().
}
  \item{Biom.tstep}{
Integer. The number of time steps over which to average the population biomass 
counting from the end of the season backwards.
}
  \item{Biom.xpos}{
Numeric. The position of the biomass value on the x-axis of the first panel of 
the plot in relative units.
}
  \item{Biom.ypos}{
Numeric. The position of the biomass value on the y-axis of the first panel of 
the plot in relative units.
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

\author{
Ruben H. Roa-Ureta
}
\note{
The target symbols on the bottom of the tope left panel are the timings of any
perturbations set by the user.
In transit fisheries, entry target symbols are in red and exit target symbols are 
in blue.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~iplot }
