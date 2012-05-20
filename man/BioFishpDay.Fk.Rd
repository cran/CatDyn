\name{BioFishpDay.Fk}
\alias{BioFishpDay.Fk}
\title{
Process Biological Sampling and Fisheries Census Data FIFD-Style Databases
}
\description{
Construction of a list containing elements necessary to fit Catch Dynamics Models
from databases maintained in the style of the Falkland Islands Fisheries Department.
}
\usage{
BioFishpDay.Fk(BioData, CatEffData, Year, Season, StartDate, EndDate,
               FGTable, span, sds, plot)
}
\arguments{
  \item{BioData}{
A data.frame produced by MBWpDaypSam().
}
  \item{CatEffData}{
A data.frame containing vessel id, date, grid square, trawl time in hours, catch in kg, fishing ground,
season, and year, for each day of fishing operation.
}
  \item{Year}{
Year of the biological and fisheries data.
}
  \item{Season}{
Either 1 (summer) or 2 (winter)
}
  \item{StartDate}{
Date of the first day of the season, in the ISO 8601 international standard yyyy-mm-dd
}
  \item{EndDate}{
Date of the last day of the season, in the ISO 8601 international standard yyyy-mm-dd.
}
  \item{FGTable}{
A table counting the number of vessel-days in the season-year in a given fishing ground.
}
  \item{span}{
A numeric value between 0 and 1 to determine the degree of smoothing in predicting missing
body mass daily data and replacing outliers using loess().
}
  \item{sds}{
Number of standard deviations away from mean to declare a mean body mass value as outlier.
}
  \item{plot}{
Logical, should a plot of daily mean body mass be returned, to see the daily values
where data was predicted because of missing data or outliers?
}
}
\details{
Use this function to create the data.frame with the input data to the models for
each season.
It calculates mean mantle length, standard deviation of mean mantle length
(yet unused), mean body mass, sd of mean body mass (yet unused), for all days
of a fishing season. It deals with missing days of biological data by prediction
using local polynomial regression (function loess). It deals with outliers by
Chebyshev inequality.
}
\value{
\item{FullSeason }{A data.frame, with data from all three fishing grounds}
\item{Beauchene }{An object of class CatDynData, a data.frame, with results for Beauchene }
\item{Central }{An object of class CatDynData, a data.frame, with results for Central }
\item{North }{An object of class CatDynData, a data.frame, with results for North }
}
\author{
Ruben H. Roa-Ureta
}
\examples{
data(SeasonDates)
data(CatEffData)
data(BioData)
StartDate.1990.S1   <- SeasonDates[1,2]
EndDate.1990.S1     <- SeasonDates[1,3]
FGTable.1990.S1     <- table(CatEffData$Fishing.ground[CatEffData$Year==1990
                             & CatEffData$Season==1
                             & CatEffData$Date <= EndDate.1990.S1
                             & CatEffData$Date >= StartDate.1990.S1])
BioData.1990.S1     <- MBWpDaypSam.Fk(BioData,Year=1990,
                                      Season=1,
                                      col1=3,
                                      col2=55,
                                      len1=4,
                                      len2=30,
                                      w1=0.14019,
                                      sd.w1=0.0017367,
                                      w2=2.2693,
                                      sd.w2=0.0050422,
                                      cor.w1.w2=-0.9876)
SeasonData.1990.S1  <- BioFishpDay.Fk(BioData=BioData.1990.S1,
                                      CatEffData,Year=1990,
                                      Season=1,
                                      StartDate=StartDate.1990.S1,
                                      EndDate=EndDate.1990.S1,
                                      FGTable=FGTable.1990.S1,
                                      span=1,
                                      sds=5,
                                      plot=TRUE)
SeasonData.1990.S1
}

\keyword{ datasets }
