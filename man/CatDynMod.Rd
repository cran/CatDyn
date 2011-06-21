\name{CatDynMod}
\alias{CatDynMod}
\title{
Generic CatDyn Results Plot
}
\description{
Generic function to plot CatDyn results objects.
}
\usage{
CatDynMod(x, ...)
}
\arguments{
  \item{x}{
An object of class CatDynMod
}
  \item{\dots}{
Not used
}
}
\details{
This generic function will just be used to give a class to the data.frame output by
CDMN0P, ..., CDMN4P.
}
\value{
A class attribute
}
\author{
Ruben H. Roa-Ureta
}
\note{
Objects of class CatDynMod are plotted by plot.CatDynMod to examine initial parameter
values before optimization and estimated parameters values after optimization.
}
\examples{
#See example in man pages for CDMN0P(), ..., CDMN4P, and plot.CatDynMod().
}