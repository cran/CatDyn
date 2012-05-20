\name{SeasonData.1990.S1.B}
\alias{SeasonData.1990.S1.B}
\docType{data}
\title{
Processed Data to Fit Catch Dynamics Models
}
\description{
Daily step, catch in billions, two measures of effort, and the mean body mass series corresponding to the summer fishing season
of the Loligo gahi fishery of the Falkland Islands
}
\usage{data(SeasonData.1990.S1.B)}
\format{
  A data frame with 120 observations on the following 5 variables.
  \describe{
    \item{\code{period}}{a numeric vector}
    \item{\code{obscat}}{a numeric vector}
    \item{\code{obseff1}}{a numeric vector}
    \item{\code{obseff2}}{a numeric vector}
    \item{\code{obsmbm}}{a numeric vector}
  }
}
\details{
See the example in BioFishpDay.Fk() to check how this database was created out of raw data from observers at sea.
}
\source{
http://fis.com/falklandfish/
}
\examples{
data(SeasonData.1990.S1.B)
head(SeasonData.1990.S1.B,10)
}
\keyword{datasets}
