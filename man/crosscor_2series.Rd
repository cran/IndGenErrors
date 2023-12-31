% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crosscor_2series.R
\name{crosscor_2series}
\alias{crosscor_2series}
\title{Cross-correlations for testing independence between the innovations of 2 series of same length}
\usage{
crosscor_2series(x, y, lag, graph = TRUE)
}
\arguments{
\item{x}{Pseudo-observations (or residuals) of first series}

\item{y}{Pseudo-observations (or residuals) of second series}

\item{lag}{Maximum number of lags around 0}

\item{graph}{Set to TRUE for a correlogram for all possible lags.}
}
\value{
\item{stat}{Cross-correlations for all lags}

\item{LB}{Sum of squares of cross-correlations}

\item{pvalue}{P-value of LB}

\item{subsets}{c(-lag:lag)}

\item{n}{length of the time series}
}
\description{
This function computes the cross-correlations  between x(t) and y(t-l), for l=-lag,.., lag, and also the combination (Wald's type)  of these statistics.
}
\examples{
data(gas)
outr <-crosscor_2series(gas$xres,gas$yres,3)

}
\references{
Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
}
