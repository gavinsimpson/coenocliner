\name{sim1dNegbinom}
\alias{sim1dNegbinom}
\title{Simulate negative binomial counts along a single gradient}
\usage{
sim1dNegbinom(x, opt, tol, h, alpha)
}
\arguments{
  \item{x}{numeric; gradient locations.}

  \item{opt}{numeric; species optima, one per taxon.}

  \item{tol}{numeric; species tolerances, one per taxon.}

  \item{h}{numeric; species abundance at optima, one per
  taxon.}

  \item{alpha}{numeric; dispersion parameter for the
  negative binomial.}
}
\value{
a matrix of simulated counts.
}
\description{
Simulates species abundances along a single gradient
simulated from a negative binomial distribution assuming a
Gaussian response model.
}
\details{
The response model used is the classical Gaussian response
with parameters \eqn{\mu}{mu}, the species optimum,
\eqn{t}, the species tolerance, and \eqn{h}, the height of
the response curve at the optimum.
}
\examples{
set.seed(1)
x1 <- runif(300, min = 4, max = 6)
Opt <- seq(4, 6, length = 5)
Tol <- rep(0.25, 5)
H <- rep(20, 5)
y <- sim1dNegbinom(x1, Opt, Tol, H, alpha = 1.1)
}
\author{
Gavin L. Simpson
}

