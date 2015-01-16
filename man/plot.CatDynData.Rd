\name{plot.CatDynData}
\alias{plot.CatDynData}
\title{
Exploratory Analysis of Catch and Effort Fisheries Data
}
\description{
Allows examining the relation between catch and effort, the marginal distributions
of catch and effort, and the time series of catch, effort, the catch spike
statistic, and mean body mass in the catch.
}
\usage{
\method{plot}{CatDynData}(x, mark, offset, hem, \dots)
}
\arguments{
  \item{x}{
An object of class CatDynData.
}
  \item{mark}{
Logical. If TRUE then the time step is posted on top of each point of the time
series of catch, effort, and the catch spike statistic.
}
  \item{offset}{
Numeric. A vector of length 3 that positions the mark above a given distance 
over the point.
}
  \item{hem}{
Character. Either "N" (northern hemisphere) or "S" (southern hemisphere).
}
  \item{\dots}{
Further arguments to be passed to plot(), hist().
}
}
\details{
Use NA to cancel the mark over the points of any of the three time series that can
be marked.
In the case of two-fleet models, the plot will display the data for the first
fleet, then the user needs to hit Enter to display the data for the second fleet.
}
\value{
A seven panel plot.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~iplot }

