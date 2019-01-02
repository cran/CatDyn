\name{CatDynBSD}
\alias{CatDynBSD}
\title{
Calculate Annual Biomass and its Standard Error
}
\description{
Using results of a model fit by CatDynFit, calculate annual biomass considering
initial abundance, mean weight, abundance inputs, fishery removals, and natural
mortality, and then computes its standard error using the delta method.
}
\usage{
CatDynBSD(x, method, multi, mbw.sd)
}
\arguments{
  \item{x}{
One object of class catdyn when the time step is the month (multi-annual version)
or a list of objects of class catdyn when the time step is the week or the day
(intra-annual versions). The number of years when the time step is the month or
the number of separate catdyn objects when the time step is the week or the day
must be 15 or higher.
}
  \item{method}{
Character or character vector. The numerical method used to fit the model when the
time step is the month or a character vector with all the methods used to fit the
several models when the time step is the week or the day.
}
  \item{multi}{
Logical. Wheter there is a single multi-annual model (monthly time step), or
many intra-annual models (daily or weekly time step)
}
  \item{mbw.sd}{
Numeric. A vector of length 12 with the standard deviation or standard error of
the mean weight in kg per month when the time step is the month, or a data.frame
with three columns when the time step is the week or the day: year, mean weight
in kg, and standard deviation or standard error of mean weight in kg. In the latter
case the number of rows must equal the number of years (>14), also the length of
the list x.
}
}
\details{
The main purpose of this function is to obtain annual biomass estimates to be passed
as input information for the fit of a population dynamics of the surplus production
kind in a hierarchical inference framework. Thus is carries over most of the
uncertainty in the original catch, effort and mean weight data, to inform the
population dynamics model. The limit of 15 years of data as a minimum for the use
of this function is set so that the fit of the population dynamics model has
sufficient information to estimate its parameters.

When the time step is the month, this function will calculate the biomass and its
standard error for every month in the time series of data. When the time step is
the week or the day it will calculate the biomass and its standard error at the
start of the season, only one value per year.
}
\value{
When the time step is the month, a data.frame with columns for year, month, time
step, mean weight (kg), standard deviation or standard error of mean weight,
abundance, standard error of abundance, biomass (tons), and standard error of
biomass (tons). The data.frame has as many rows as time steps (minimum of 180 months)

When the time step is the week or the day, a data.frame with the year, mean weight
(kg), standard deviation or standard error of mean weight, abundance, standard error
of abundance, biomass (tons), and standard error of biomass (tons). The data.frame
has as many rows as years (minimum of 15 years).
}
\references{
Roa-Ureta, R. H. et al. 2015. Fisheries Research 171 (Special Issue), 59-67.
}
\author{
Ruben H. Roa-Ureta (ORCID ID 0000-0002-9620-5224)
}
\note{
This function makes extensive use of the delta method for carrying over the original
uncertainty. It uses the asymptotic standard errors of parameters in the object of
class catdyn and their correlation matrix, along with standard deviation or
standard error of mean weight per time step, to calculate the standard error of
biomass.

In case a model output (object of class catdyn) is selected which has not produced
standard errors for all parameters that are involved in the calculation of biomass
(natural mortality, initial abundance, and magnitude of each input/output pulses of
abundance) then the missing standard errors are replaced by imputed standard errors
computed as the estimate for which the standard error is missing times the mean
coefficient of variation across all parameters which did get standard errors. In this
manner the degree of statistical uncertainty is preserved and all standard errors
are available to use the delta method to calculate the standard error of biomass.
}
\keyword{ ~manip }
\keyword{ ~optimize }
\keyword{ ~models }

