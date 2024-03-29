% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conserved.R
\name{mirPlot}
\alias{mirPlot}
\title{mirPlot}
\usage{
mirPlot(dfm, xTxtSize = 10, yTxtSize = 10, plotWidth = 170, plotHeight = 150)
}
\arguments{
\item{dfm}{Input dataframe.}

\item{xTxtSize}{Font size of x-axis tick label. Default is 10.}

\item{yTxtSize}{Font size of y-axis tick label. Default is 10.}

\item{plotWidth}{Width of the output image file. Default is 170.}

\item{plotHeight}{Height of the output image file. Default is 150.}
}
\value{
Histogram of the miRNA representative expression levels
}
\description{
Histogram. Also usable in most of the situations, expecially useful when end = 1.
}
\examples{
\dontrun{
mirPlot(miRNANrm, xTxtSize = 6, plotWidth = 360)
}
}
