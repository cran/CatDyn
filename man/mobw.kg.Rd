\name{mobw.kg}
\alias{mobw.kg}
\title{
Mean body weight per month
}
\description{
Using a data.frame with paired observations of length and weight of individuals, and
a second data.frame with paired observations of length and month of individuals, predict
the mean weight per month in kg and its standard error.
}
\usage{
mobw.kg(par, lenw, lenm, method, span)
}
\arguments{
  \item{par}{
Numeric vector length 2 with the initial values for parameters of the length-weight relationship.
}
  \item{lenw}{
Numeric data.frame with paired observations of length and weight (kg).
}
  \item{lenm}{
Numeric data.frame with paired observations of length and month in numeric format.
}
  \item{method}{
Character, the numerical method to fit the length-weight relationship. Use one of
the following: BFGS,CG,Nelder-Mead,SANN.
}
  \item{span}{
The span parameter of the loess function to fit the length-month relationship.
}
}
\details{
The length data in lenw and lenm must be in the same units.
}
\value{
A data.frame with five columns: month, predicted mean weight from the length data
and the length-weight relationship, its standard error, smoothed mean length from
the loess fit of the length-month relationship, and its standard error.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\examples{
data(gayhakelm)
data(gayhakelw)
par <-  c(4e-4,2.9)
method <- "SANN"
span   <- 1
chihake.mobw <- mobw.kg(par=par, lenw=gayhakelw, lenm=gayhakelm,
                        method=method, span=span)
}
\keyword{ ~models }