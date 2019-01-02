\name{CatDynCor}
\alias{CatDynCor}

\title{
Correlation Plot Among Parameter Estimates
}
\description{
It extracts the correlation matrices output by CatDynFit() from several model fits
and plots histograms of pairwise correlation coefficients.
}
\usage{
CatDynCor(x, ttl, method, arr)
}
\arguments{
  \item{x}{
A list of model objects of class catdyn output by CatDynFit().
}
  \item{ttl}{
A character vector to use as title of each histogram.
}
  \item{method}{
A character vector of numerical methods to extract correlation matrices.
}
  \item{arr}{
A numerical vector of length 2 used to organize the panels containing each histogram.
Passed as is to par().
}
}
\details{
The arguments x, ttl, and method must be of the same length.

It might be useful to examine the results of different numerical methods applied
to the same model with this function. To do this just repeat the name of the 
object in x and specify the different numerical methods in method.
}
\value{
A multiple panel correlation plot.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
Those histograms that show more correlations concentrated around zero indicate a 
better fit.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~graphics }
\keyword{ ~models }
