\name{biGaussianResponse}
\alias{biGaussianResponse}
\title{Bivariate Gaussian response model}
\usage{
biGaussianResponse(x1, x2, opt1, tol1, opt2, tol2, h, corr = 0)
}
\arguments{
  \item{x1}{numeric; locations along gradient 1.}

  \item{x2}{numeric; locations along gradient 2.}

  \item{opt1}{numeric; species optima on gradient 1, one
  per taxon.}

  \item{tol1}{numeric; species toleranceon gradient 1, one
  per taxon.}

  \item{h}{numeric; species abundance at the optimum}

  \item{opt2}{numeric; species optima on gradient 2, one
  per taxon.}

  \item{tol2}{numeric; species tolerance on gradient 2, one
  per taxon.}

  \item{corr}{numeric; the correlation between gradients.}
}
\value{
A vector of expected responses according to the bivariate
Gaussian repsonse model.
}
\description{
Bivariate Gaussian response model
}
\author{
Gavin Simpson
}

