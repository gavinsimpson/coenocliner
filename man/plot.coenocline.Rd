% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plot.coenocline.R
\name{plot.coenocline}
\alias{lines.coenocline}
\alias{plot.coenocline}
\title{Plot species simulations along gradients}
\usage{
\method{plot}{coenocline}(x, type = "p", ...)

\method{lines}{coenocline}(x, ...)
}
\arguments{
\item{x}{an object of class \code{"coenocline"}, the result of a call to \code{\link{coenocline}}.}

\item{type}{character; the type of plot to produce. See \code{\link{plot.default}} for details.}

\item{...}{additional arguments to \code{\link{matplot}}.}
}
\value{
A plot is drawn on the current device.
}
\description{
A simple S3 \code{\link{plot}} method for coenocline simulations.
}
\author{
Gavin L. Simpson
}
\keyword{hplot}
