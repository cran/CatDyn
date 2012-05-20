\name{MBWpDaypSam.Fk}
\alias{MBWpDaypSam.Fk}
\title{
Mean Bodymass per Day and Sample
}
\description{
This function estimates the MLEs of mean body mass per day and their standard deviations
from the length frequencies and the estimates of parameters of a length-bodymass power relationship,
including the variances and covariances of these estimates. Custom-made for FIFD databases.
}
\usage{
MBWpDaypSam.Fk(BioData, Year, Season, col1, col2, len1, len2, w1,
               sd.w1, w2, sd.w2, cor.w1.w2)
}
\arguments{
  \item{BioData}{
A data.frame consisting at least of a year column, a season column, and several columns containing length frequencies.
}
  \item{Year}{
Year of the biological and fisheries data.
}
  \item{Season}{
Either 1 (summer) or 2 (winter)
}
  \item{col1}{
Column index of the first length class, 14 cm.
}
  \item{col2}{
Column index of the last length class, 30 cm.
}
  \item{len1}{
The smaller length class, 14.
}
  \item{len2}{
The largest length class, 30.
}
  \item{w1}{
MLE of the linear parameter of the length-weight relationship.
}
  \item{sd.w1}{
Estimated standard deviation of MLE of the linear parameter of the length-weight relationship.
}
  \item{w2}{
MLE of the power parameter of the length-weight relationship.
}
  \item{sd.w2}{
Estimated standard deviation of MLE of the power parameter of the length-weight relationship.
}
  \item{cor.w1.w2}{
Estimated covariance matrix of MLEs of the length-weight relationship.
}
}
\details{
The algorithm for standard deviations is based on a Taylor series approximation.
}
\value{
A BioData data.frame augmented with four new columns: mean length per day, standard deviation of mean
length per day, mean body mass per day, standard deviation of mean body mass per day.
}
\author{
Ruben H. Roa-Ureta
}
\examples{
#Data - 1st Season, Beauchene area, 1990, Loligo gahi fishery in the Falklands
data(SeasonDates)
data(CatEffData)
data(BioData)
StartDate.1990.S1   <- SeasonDates[1,2]
EndDate.1990.S1     <- SeasonDates[1,3]
FGTable.1990.S1     <- table(CatEffData$Fishing.ground[
                              CatEffData$Year==1990 &
                              CatEffData$Season==1 &
                              CatEffData$Date <= EndDate.1990.S1 &
                              CatEffData$Date >= StartDate.1990.S1])
#Next process can take a minute or two ...
MBWpDaypSam.Fk(BioData,
               Year=1990,
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
}
\keyword{ datasets }
