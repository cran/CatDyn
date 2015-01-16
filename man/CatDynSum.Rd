\name{CatDynSum}
\alias{CatDynSum}
\title{
Summary Table for Model Comparison
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
A character vector of numerical methods to extract from each element of the 'x' list.
}
}
\details{
The arguments 'x' and 'meths' must be of the same length.
If a single model fit object of class 'catdyn' was fit with more than one numerical method and you
wish to include a row of the summary for more than one of those methods, then repeat the name of the
object in 'x' and specify the numerical methods in 'meths'.
}
\value{
A data.frame with an extraction of information from several models fitted to the data.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~misc }
