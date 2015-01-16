\name{CatDynCor}
\alias{CatDynCor}

\title{
Correlation Plot Among Parameter Estimates
}
\description{
It extracts the correlation matrices output by optimx() numerical methods from several model fits
and plots histograms of pairwise correlation coefficients.
}
\usage{
CatDynCor(x, ttl, meths, arr)
}
\arguments{
  \item{x}{
A list of model objects of class 'catdyn' output by CatDynFit().
}
  \item{ttl}{
A character vector to use as title of each histogram.
}
  \item{meths}{
A character vector of numerical methods to extract correlation matrices.
}
  \item{arr}{
A numerical vector of length 2 used to organize the panels containing each histogram. Passed as is to par().
}
}
\details{
The arguments 'x', 'ttl', and 'meths' must be of the same length.
If a single model fit object of class 'catdyn' was fit with more than one numerical method and you
wish to plot the correlation histogram plot for more than one of those methods, then repeat the name of the
object in 'x' and specify the numerical methods in 'meths'.
}
\value{
A multiple panel correlation plot.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~graphics }
