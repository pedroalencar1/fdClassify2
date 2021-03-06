% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpentad.R
\name{f.pentad}
\alias{f.pentad}
\title{Accumulation into pentads}
\usage{
f.pentad(vtime, vvalue, na_rm = F, f = mean)
}
\arguments{
\item{vtime}{data frame column or vector containing \code{date} data}

\item{vvalue}{data frame column or vector containing the analysed data}

\item{na_rm}{boolean (should NA values be removed? Defaulf = F)}

\item{f}{\code{R function} to be applied (default = mean)}
}
\value{
The function return a list with two elements. One data frame with time stamped pentad values and a matrix with the 73 pentads organized in lines and years in columns.
}
\description{
Internal function to accumulate data into pentads (5-day long periods) using different accumulation functions (mean, max, min, sum, etc.)
}
