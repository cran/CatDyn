\name{plot.CatDynData}
\alias{plot.CatDynData}
\title{
Exploratory analysis of catch and effort fisheries data on a graphical display
}
\description{
A six panel plot of raw catch, effort, and biological sampling per time step.
}
\usage{
plot.CatDynData(x, tstep, mult, unit1, unit2, bmunit, span,
                top.text, hem, ...)
}
\arguments{
  \item{x}{
A data.frame, an object of class CatDynData.
}
  \item{tstep}{
Character, days recommended.
}
  \item{mult}{
Character, the abundance multiplier.
}
  \item{unit1}{
Character, the unit of one measure of effort.
}
  \item{unit2}{
Character, the unit of a second measure of effort.
}
  \item{bmunit}{
Character, the unit of body mass.
}
  \item{span}{
A numeric value between 0 and 1 to determine the degree of smoothing in predicting missing
body mass daily data and replacing outliers using loess().
}
  \item{top.text}{
Character, a descriptive sentence.
}
  \item{hem}{
Character, either N for northern hemisphere or S for southern hemisphere.
}
  \item{\dots}{
Arguments to be passed to methods.
}
}
\details{
The hem argument allows writing the season of the year over the plots.

The plot is mostly useful to observe which of two effort measures better relates
to catch, but it also shows along-season trends of catch, effort, and mean body
mass, the latter useful to check for the presence of perturbations due to
immigration of small, juvenile fish.
}
\value{
A six-panel plot.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonData.1990.S1.B)
class(SeasonData.1990.S1.B) <- "CatDynData"
#Plot
plot(x=SeasonData.1990.S1.B,
     tstep="Day",
     mult='Billions',
     unit1="No Vessels",
     unit2="Hours",
     bmunit="kg",
     span=1,
     top.text="Beauchene Season 1 1990",
     hem='S')
}
\keyword{ optimize }
\keyword{ models }
\keyword{ & &nonlinear  }
