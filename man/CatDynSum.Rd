\name{CatDynSum}
\alias{CatDynSum}
\title{
Summary Table for Generalized Depletion (Catch Dynamics) Model Comparison
}
\description{
Summarize information theoretic, numerical and statistical results of several models
fit to the data to select a best working model.
}
\usage{
CatDynSum(x, season, meths)
}
\arguments{
  \item{x}{
A list of model objects of class 'catdyn' output by CatDynFit().
}
  \item{season}{
Character, an identifier for the set of models.
}
  \item{meths}{
A character vector of numerical methods to extract from each element of the 'x' 
list.
}
}
\details{
The arguments 'x' and 'meths' must be of the same length.

It might be useful to examine the results of different numerical methods applied
to the same model with this function. To do this just repeat the name of the 
object in 'x' and specify the different numerical methods in 'meths'.
}
\value{
A data.frame with an extraction of comparative information from several models 
fitted to the data.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~misc }
