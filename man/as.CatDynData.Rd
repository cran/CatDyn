\name{as.CatDynData}
\alias{as.CatDynData}
\title{
Data Object for the Estimation of Catch Dynamic Models
}
\description{
It takes the vectors of catch, effort, and mean body weight from a dataframe and 
creates an object of class CatDynData. Objects of this class are lists with two 
components, one for properties of the data such as units and another for the data: 
catch, effort, mean body weight by fleet, and the catch spike statistic.
}
\usage{
as.CatDynData(x, step, fleet.name, coleff, colcat, colmbw, 
              unitseff, unitscat, unitsmbw, nmult, season.dates)
}
\arguments{
  \item{x}{
A dataframe where to find the columns of catch, effort, and mean body weight
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
Integer. The column(s) in "x" where to find the mean body weight data.
}
  \item{unitseff}{
Character. The unit(s) of effort.
}
  \item{unitscat}{
Character. The unit of catch. Either "ton" (metric tonnes), "kg", or "ind" 
(individuals).
}
  \item{unitsmbw}{
Character. The unit of body weight. Either "kg", "g", or "ind" (individuals).
}
  \item{nmult}{
Character. The multiplier that scales the catch in numbers. Either "bill" 
(billions), "mill" (millions), "thou" (thousands), or "ind" (individuals).
}
  \item{season.dates}{
Character vector. A two component character vector with the initial and final 
dates of the season in the ISO 8601 standard.
}
}
\details{
The time step determines the rows of 'x'. Make it sure that the number of rows, 
i.e. the length of the season in time steps, is large enough to estimate all 
parameters in the model. The simplest model has five parameters, the most complex 
model has 50 parameters. A rule of thumb would be that the number of time steps 
be at least three times the number of parameters. 

If it is a two fleet system, combine the fleet names, such as c("industrial", 
"artisanal"), and likewise with coleff, colcat, and colmbw, such as c(5,9) to 
indicate the columns of catch for the industrial and artisanal fleets respectively. 
The same applies to units of effort. In a two fleet system, the time step, and 
the units of catch, mean body weight, and the multiplier must be the same for both 
fleets. 

When the unit of catch and of body weight is "ind", it means that the catch was 
counted in numbers, not in biomass. In that case the mandatory column of mean body 
weight should be a column of 1s. The multiplier is the quantity by which the catch 
in the model shall be raised to be scaled to the actual catch. The idea here is 
that in many fisheries the daily, weekly, or monthly catch (for example, anchovies, 
squids) is very large so by setting the multiplier to "bill", "mill", or "thou", 
the model is working with catches in the orders of tens at most. If the multiplier 
is set to "ind" then the catch is modeled at the level of the actual catch by time 
step. This option is useful for sport fisheries, in combination with the 'poisson' 
or 'negbin' option for distribution.

The 'season.dates' parameter will allow counting the number of steps in integer 
sequential fashion. If the "time.step" parameter is "day" or "week" the dates may 
jump one year at most, whereas if the 'time.step' parameter is "month" then 
'season.dates' may jump over many years. When the time step is week or months, 
this parameter needs not be precisely specified; any day within the right week 
or month will suffice. If you get an error message saying that the number of time 
steps is not right and that you should consider changing 'season.dates', then just 
change the dates a few days until you no longer get the error message.

The catch spike statistic is a fleet-specific statistic that is useful to identify 
the timing of perturbations to depletion; it is defined as
\deqn{S_{f,t} = 10 \times (X_{f,t}/max(X_{f,t}) - E_{f,t}/max(E_{f,t}))}
where \emph{X} is the observed catch at time step \emph{t} by fleet \emph{f}, and 
\emph{E} is the observed effort. When this statistic is positive and high then the 
time step at which this happened is a candidate for a perturbation. In transit 
fisheries the complementary reasoning is valid: when the statistic is negative and 
low then the time step at which this happened is a candidate for the timing of one 
emmigration event.
}
\value{
A list of length 2.
\item{Properties }{The units for time step, catch, body weight, and the catch 
                   numbers multiplier; the names of the fleets and the units of 
                   effort, and the start date and end date of the fishing season}
\item{Data }{One dataframe for each fleet with the time step, effort, catch in 
             biomass, mean body weight, catch in the numbers multiplier, and the
             catch spike statistic}

}
\references{
Roa-Ureta, R. H. 2012. ICES J. Mar. Sci. 69:1403-1415.
Roa-Ureta, R. H. 2015. Fisheries Research (In Press), doi:10.1016/j.fishres.2014.08.014
Roa-Ureta, R. H. et al. 2015. Fisheries Research (In Press), doi:10.1016/j.fishres.2014.12.006
}
\author{
Ruben H. Roa-Ureta
}
\note{
The objects created with as.CatDynData will pass the raw data to plotting and 
estimating functions.
}
\examples{
lgahi <- as.CatDynData(x=lolgahi,
                       step="day",
                       fleet.name="Fk",
                       coleff=2,
                       colcat=1,
                       colmbw=3,
                       unitseff="nboats",
                       unitscat="kg",
                       unitsmbw="kg",
                       nmult="bill",
                       season.dates=c(as.Date("1990-01-31"),
                                      as.Date("1990-05-30")))
}
\keyword{ ~manip }
