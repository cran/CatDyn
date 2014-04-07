\name{deltamethod}
\alias{deltamethod}
\title{
First-order Taylor Series Expansion of Functions of Random Variables
}
\description{
The delta method for approximating the standard error of a transformation \emph{g}(\emph{X}) 
of a random variable \emph{X} = (\emph{x}1, \emph{x}2, ...), given estimates of the mean and covariance 
matrix of \emph{X}. 
}
\usage{
deltamethod(g, mean, cov, ses = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
A formula representing the transformation. The variables must be 
labelled \emph{x}1, \emph{x}2, ... For example, ~ 1 / (\emph{x}1 + \emph{x}2). If the transformation returns 
a vector, then a list of formulae representing (\emph{g}1, \emph{g}2, ...) can be provided, for 
example list( ~ \emph{x}1 + \emph{x}2, ~ \emph{x}1 / (\emph{x}1 + \emph{x}2) ).
}
  \item{mean}{
The estimated mean of \emph{X}.
}
  \item{cov}{
The estimated covariance matrix of \emph{X}.
}
  \item{ses}{
If TRUE, then the standard errors of \emph{g}1(\emph{X}), \emph{g}2(\emph{X}), ... are returned. Otherwise the 
covariance matrix of \emph{g}(\emph{X}) is returned.
}
}
\details{
This function was copied from package msm. It is used in CatDyn to backtransform
from the logarithm because CatDyn parameters are all estimated in the log scale
to improve numerical performance.
For more details see the help pages for function deltamethod of package msm.
}
\value{
A vector containing the standard errors of \emph{g}1(\emph{X}), \emph{g}2(\emph{X}), ... or a matrix containing the covariance of \emph{g}(\emph{X}). 
}
\author{
C. H. Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\examples{
#See the examples in package msm.
}