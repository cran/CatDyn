\name{M.Hoenig}
\alias{M.Hoenig}
\title{
Estimate natural mortality rate from longevity data
}
\description{
Using the empirical regression parameters obtained by Prof. John M. Hoenig
and the longevity provided with max.age, estimate the natural mortality rate and its
standard error from application of the delta method.
}
\usage{
M.Hoenig(max.age, time.step)
}
\arguments{
  \item{max.age}{
Numeric, longevity in years.
}
  \item{time.step}{
Character, either day, week, or month, the desired unit of the estimate.
}
}
\details{
This function can be used to obtain a good initial value for natural mortality
rate to be input in the exploratory or statistical feet of a catch dynamics model,
at the time step in which the model will be run.
}
\value{
A length 2 numeric vector with estimated natural mortality and its standard error.
}
\references{
Hoening, J.M. 2005. Empirical use of longevity data to estimate mortality rates.
SEDAR33-RD17. SEDAR, North Charleston, SC. 8 pp.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\examples{
max.age <- 5.8
time.step <- "day"
M.Hoenig(max.age,time.step)
}
\keyword{ ~models }
