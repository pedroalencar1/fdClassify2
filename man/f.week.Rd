% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fweek.R
\name{f.week}
\alias{f.week}
\title{Accumulation into weeks}
\usage{
f.week(data.var, na_rm = F, f = mean, kind = "standard")
}
\arguments{
\item{data.var}{data frame with time stamps}

\item{na_rm}{boolean (should NA values be removed? Defaulf = F)}

\item{f}{\code{R function} to be applied (default = mean)}

\item{kind}{String indicating what is the size of the week. The default value is "standard" (conventional 7 day week). Other option is "noguera", that devides each month in 4 weeks (1 to 8; 9 to 15; 16 to 22; and 23 to 29/30/31)}
}
\value{
The function return a list with two elements. One data frame with time stamped pentad values and a matrix with the 52 (or 48) weeks organized in lines and years in columns.
}
\description{
Internal function to accumulate data into weeks using different accumulation functions (mean, max, min, sum, etc.)
}
