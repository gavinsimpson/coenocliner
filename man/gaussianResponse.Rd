\name{gaussianResponse}
\alias{gaussianResponse}
\title{Univariate Gaussian response model}
\usage{
gaussianResponse(x, opt, tol, h)
}
\arguments{
  \item{x}{numeric; gradient locations.}

  \item{opt}{numeric; species optima, one per taxon.}

  \item{tol}{numeric; species tolerances, one per taxon.}

  \item{h}{numeric; species abundance at optima, one per
  taxon.}
}
\value{
A vector of expected responses according to the univariate
Gaussian response model.
}
\description{
Univariate Gaussian response model
}
\author{
Gavin Simpson
}

