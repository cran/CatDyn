\name{CatEffData}
\alias{CatEffData}
\docType{data}
\title{
Catch and Effort Census Data
}
\description{
Exhaustive recorded operational activity of the Loligo fleet off the Falkland Islands
during the first season of 1990.
}
\usage{data(CatEffData)}
\format{
  A data frame with 3774 observations on the following 8 variables.
  \describe{
    \item{\code{Vessel}}{a factor with levels \code{VV01} \code{VV02} \code{VV03} \code{VV04} \code{VV05} \code{VV06} \code{VV07} \code{VV08} \code{VV09} \code{VV10} \code{VV11} \code{VV12} \code{VV13} \code{VV14} \code{VV15} \code{VV16} \code{VV17} \code{VV18} \code{VV19} \code{VV20} \code{VV21} \code{VV22} \code{VV23} \code{VV24} \code{VV25} \code{VV26} \code{VV27} \code{VV28} \code{VV29} \code{VV30} \code{VV31} \code{VV32} \code{VV33} \code{VV34} \code{VV35} \code{VV36} \code{VV37} \code{VV38} \code{VV39} \code{VV40} \code{VV41} \code{VV42} \code{VV43} \code{VV44} \code{VV45} \code{VV46} \code{VV47} \code{VV48} \code{VV49}}
    \item{\code{Date}}{a Date}
    \item{\code{Square}}{a factor with levels \code{XLAN} \code{XLAP} \code{XMAP} \code{XMAQ} \code{XNAP} \code{XNAQ} \code{XPAP} \code{XQAN} \code{XQAP} \code{XRAN} \code{XRAP} \code{XSAM} \code{XSAN} \code{XSAP} \code{XTAG} \code{XTAH} \code{XTAL} \code{XTAM} \code{XTAN} \code{XUAG} \code{XUAH} \code{XUAJ} \code{XUAK} \code{XUAL} \code{XUAM} \code{XUAN} \code{XVAG} \code{XVAH} \code{XVAJ} \code{XVAK} \code{XVAL} \code{XVAM} \code{XVAN} \code{XVAP} \code{XVAQ} \code{XWAH} \code{XWAJ} \code{XWAK} \code{XWAL}}
    \item{\code{Trawl.Time.h}}{a numeric vector}
    \item{\code{Catch.kg}}{a numeric vector}
    \item{\code{Fishing.ground}}{a factor with levels \code{Beauchene} \code{Central} \code{North}}
    \item{\code{Season}}{a numeric vector}
    \item{\code{Year}}{a numeric vector}
  }
}
\source{
http://fis.com/falklandfish/
}
\examples{
data(CatEffData)
head(CatEffData, 10)
}
\keyword{datasets}
