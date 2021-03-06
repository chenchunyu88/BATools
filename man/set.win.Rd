% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/window.R
\name{set.win}
\alias{set.win}
\title{create fix size window}
\usage{
set.win(map, len = NULL, unit = c("Mb", "kb", "count"))
}
\arguments{
\item{map}{genomic map has at least two column with column name \code{chr} for chromosome id and \code{pos} for postition in the chromosome in base pair (bp)}

\item{len}{window length}

\item{unit}{the unit for window length can be either Mb or kb or count}
}
\value{
a map with new column \code{idw}
}
\description{
create fix size window
}
\examples{
\dontrun{
 rm(list=ls())
 library(BATools)
 data("Pig")
 
}
}
