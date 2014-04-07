\name{as.CatDynData}
\alias{as.CatDynData}
\title{
Data Object for the Estimation of Catch Dynamic Models
}
\description{
It takes the vectors of catch, effort, and mean body mass from a dataframe and creates an object of class CatDynData.
Objects of class CatDynData are lists with two components, one for properties of the data such as units and another for the data: 
catch, effort, mean body mass by fleet, and the catch spike statistic.
}
\usage{
as.CatDynData(x, step, fleet.name, coleff, colcat, colmbw, 
              unitseff, unitscat, unitsmbw, nmult, season.dates)
}
\arguments{
  \item{x}{
A dataframe where to find the columns of catch, effort, and mean body mass
}
  \item{step}{
Character. The time step of the dynamics, either "day", "week", or "month".
}
  \item{fleet.name}{
Character. The name of the fleet(s). 
}
  \item{coleff}{
Integer. The column(s) in "x" where to find the effort data. 
}
  \item{colcat}{
Integer. The column(s) in "x" where to find the catch data. 
}
  \item{colmbw}{
Integer. The column(s) in "x" where to find the mean body mass data.
}
  \item{unitseff}{
Character. The unit(s) of effort.
}
  \item{unitscat}{
Character. The unit of catch. Either "ton" (metric tonnes), "kg", or "ind" (individuals).
}
  \item{unitsmbw}{
Character. The unit of body mass. Either "kg", "g", or "ind" (individuals).
}
  \item{nmult}{
Character. The multiplier that scales the catch in numbers. Either "bill" (billions), "mill" (millions), "thou" (thousands), or "ind" (individuals).
}
  \item{season.dates}{
Character vector. A two component character vector with the initial and final dates of the season in the ISO 8601 standard.
}
}
\details{
The time step determines the rows of 'x'. Make it sure that the number of rows, i.e. the length of the season in time steps, 
is large enough to estimate all parameters in the model. The largest model may have 48 parameters so it would be wise 
to have a few hundred time steps.
If it is a two fleet system, combine the fleet names, such as c("industrial", "artisanal"), and likewise with coleff, colcat, and colmbw,
such as c(5,9) to indicate the columns of catch for the industrial and artisanal fleets respectively. The same applies to units of effort.
In a two fleet system, the time step, and the units of catch, mean body mass, and the multiplier must be the same for both fleets.
When the unit of catch and of body mass is "ind", it means that the catch was counted in numbers, not in biomass. In that case the mandatory 
column of mean body mass should be a column of 1s.
The multiplier is the quantity by which the catch in the model shall be raised to be scaled to the actual catch. The idea here is that
in many fisheries the daily, weekly, or monthly catch (for example, anchovies, squids) is very large so by setting the multiplier to "bill" the
model is working with catches in the orders of tens at most. If the multiplier is set to "ind" then the catch is modeled at the
level of the actual catch by time step. This option is useful for sport fisheries.
The 'season.dates' parameter will allow counting the number of steps in integer sequential fashion. If the "time.step" parameter 
is "day" or "week" the dates may jump one year at most, whereas if the 'time.step' parameter is "month" then 'season.dates' may 
jump over many years.
If the time step is month, use the first day of the first month as start date, and the last day of the last month as end date
}
\value{
A list of length 2.
\item{Properties }{The units for time step, catch, body mass, and the catch numbers multiplier; the names of the fleets and the units of effort,
                   and the start date and end date of the fishing season}
\item{Data }{One dataframe for each fleet with the time step, effort, catch in biomass, mean body mass, catch in the numbers multiplier, and the
             catch spike statistic}
The catch spike statistic is a fleet-specific statistic that is useful to identify the timing of perturbations to depletion; it is defined as
\deqn{S_{f,t} = 10 \times (X_{f,t}/max(X_{f,t}) - E_{f,t}/max(E_{f,t}))}
where \emph{X} is the observed catch at time step \emph{t} by fleet \emph{f}, and \emph{E} is the observed effort. 
}
\references{
Roa-Ureta, R. H. 2012. ICES J. Mar. Sci. 69:1403-1415.
}
\author{
Ruben H. Roa-Ureta
}
\note{
The objects created with as.CatDynData will pass the raw data to plotting and estimating functions.
}
\examples{
#See examples for CatDynFit().
}
\keyword{ ~manip }
