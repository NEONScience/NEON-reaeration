% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{reaFormatted}
\alias{reaFormatted}
\title{Formatted reaeration data.}
\format{
A data frame with 24 rows and 16 variables:
\describe{
  \item{siteID}{site of water sample collection, string}
  \item{namedLocation}{location of water sample collection, string}
  \item{startDate}{date and time of experiment, dateTime}
  \item{stationToInjectionDistance}{distance between injection location and sampling location, meter}
  \item{injectionType}{either constant rate or slug, string}
  \item{slugPourTime}{date and time of slug injection, dateTime}
  \item{dripStartTime}{start date and time of constant rate injection, dateTime}
  \item{backgroundSaltConc}{background stream salt concentration, milligramsPerLiter}
  \item{plateauSaltConc}{salt concentration at plateau, milligramsPerLiter}
  \item{corrPlatSaltConc}{plateau salt concentration corrected for background salt concentration, milligramsPerLiter}
  \item{plateauGasConc}{gas concentration at plateau, ppmv}
  \item{wettedWidth}{mean wetted width based off of up to 30 measurements, meter}
  \item{waterTemp}{water temperature, celsius}
  \item{hoboSampleID}{identifier for logger file, string}
  \item{discharge}{stream discharge, literPerSecond}
  \item{eventID}{identifier for tracer injection experiment, string}
  ...
}
}
\source{
\url{http://data.neoninc.org/home/}
}
\usage{
reaFormatted
}
\description{
A dataset containing example data for calculating reaeration rates.
}
\keyword{datasets}
