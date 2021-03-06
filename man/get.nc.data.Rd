% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getncdata.R
\name{get.nc.data}
\alias{get.nc.data}
\title{get.nc.data - Function to extract data from ERA5 raw data}
\usage{
get.nc.data(my_lon, my_lat, my_filename, vname, file = T)
}
\arguments{
\item{my_lon}{Longituge of study area (decimal degrees)}

\item{my_lat}{Latitude of study area (decimal degrees)}

\item{my_filename}{File name}

\item{vname}{Variable to be extracted}

\item{file}{boolean, indicate if a csv file should be generated (defaulf = True)}
}
\value{
The function returns a data frame containing one variable and time stamps
}
\description{
get.nc.data - Function to extract data from ERA5 raw data
}
