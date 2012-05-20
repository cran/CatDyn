\name{BioData}
\alias{BioData}
\docType{data}
\title{
Biological sampling data
}
\description{
Length frequency data collected by observers onboard by day of fishing operation
and by spatial grid in the Loligo gahi fishery off the east coast of the
Falkland Islands, in half cm bins from 14 to 30 cm of mantle length, year 1990, summer season
}
\usage{data(BioData)}
\format{
  A data frame with 486 observations on the following 62 variables.
  \describe{
    \item{\code{Vessel}}{a factor with levels \code{VV01} \code{VV02} \code{VV03} \code{VV04} \code{VV05} \code{VV06} \code{VV07} \code{VV08} \code{VV09} \code{VV10} \code{VV11} \code{VV12}}
    \item{\code{Date}}{a Date}
    \item{\code{ml4cm}}{a numeric vector}
    \item{\code{ml4.5cm}}{a numeric vector}
    \item{\code{ml5cm}}{a numeric vector}
    \item{\code{ml5.5cm}}{a numeric vector}
    \item{\code{ml6cm}}{a numeric vector}
    \item{\code{ml6.5cm}}{a numeric vector}
    \item{\code{ml7cm}}{a numeric vector}
    \item{\code{ml7.5cm}}{a numeric vector}
    \item{\code{ml8cm}}{a numeric vector}
    \item{\code{ml8.5cm}}{a numeric vector}
    \item{\code{ml9cm}}{a numeric vector}
    \item{\code{ml9.5cm}}{a numeric vector}
    \item{\code{ml10cm}}{a numeric vector}
    \item{\code{ml10.5cm}}{a numeric vector}
    \item{\code{ml11cm}}{a numeric vector}
    \item{\code{ml11.5cm}}{a numeric vector}
    \item{\code{ml12cm}}{a numeric vector}
    \item{\code{ml12.5cm}}{a numeric vector}
    \item{\code{ml13cm}}{a numeric vector}
    \item{\code{ml13.5cm}}{a numeric vector}
    \item{\code{ml14cm}}{a numeric vector}
    \item{\code{ml14.5cm}}{a numeric vector}
    \item{\code{ml15cm}}{a numeric vector}
    \item{\code{ml15.5cm}}{a numeric vector}
    \item{\code{ml16cm}}{a numeric vector}
    \item{\code{ml16.5cm}}{a numeric vector}
    \item{\code{ml17cm}}{a numeric vector}
    \item{\code{ml17.5cm}}{a numeric vector}
    \item{\code{ml18cm}}{a numeric vector}
    \item{\code{ml18.5cm}}{a numeric vector}
    \item{\code{ml19cm}}{a numeric vector}
    \item{\code{ml19.5cm}}{a numeric vector}
    \item{\code{ml20cm}}{a numeric vector}
    \item{\code{ml20.5cm}}{a numeric vector}
    \item{\code{ml21cm}}{a numeric vector}
    \item{\code{ml21.5cm}}{a numeric vector}
    \item{\code{ml22cm}}{a numeric vector}
    \item{\code{ml22.5cm}}{a numeric vector}
    \item{\code{ml23cm}}{a numeric vector}
    \item{\code{ml23.5cm}}{a numeric vector}
    \item{\code{ml24cm}}{a numeric vector}
    \item{\code{ml24.5cm}}{a numeric vector}
    \item{\code{ml25cm}}{a numeric vector}
    \item{\code{ml25.5cm}}{a numeric vector}
    \item{\code{ml26cm}}{a numeric vector}
    \item{\code{ml26.5cm}}{a numeric vector}
    \item{\code{ml27cm}}{a numeric vector}
    \item{\code{ml27.5cm}}{a numeric vector}
    \item{\code{ml28cm}}{a numeric vector}
    \item{\code{ml28.5cm}}{a numeric vector}
    \item{\code{ml29cm}}{a numeric vector}
    \item{\code{ml29.5cm}}{a numeric vector}
    \item{\code{ml30cm}}{a numeric vector}
    \item{\code{Fishing.Ground}}{a factor with levels \code{Beauchene} \code{Central} \code{North}}
    \item{\code{Season}}{a numeric vector}
    \item{\code{Year}}{a numeric vector}
    \item{\code{Mlen.cm}}{a numeric vector}
    \item{\code{Slen.cm}}{a numeric vector}
    \item{\code{Mbm.kg}}{a numeric vector}
    \item{\code{Sbm.kg}}{a numeric vector}
  }
}
\details{
Usually one or two at sea observers sending information to Fisheries HQ with
a few days missing due to observer swaps.
}
\source{
http://fis.com/falklandfish/
}
\examples{
data(BioData)
head(BioData, 10)
}
\keyword{datasets}
