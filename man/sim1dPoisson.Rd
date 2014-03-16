\name{sim1dPoisson}
\alias{sim1dPoisson}
\title{Simulate Poisson counts along a single gradient}
\usage{
sim1dPoisson(x, opt, tol, h)
}
\arguments{
  \item{x}{numeric; gradient locations.}

  \item{opt}{numeric; species optima, one per taxon.}

  \item{tol}{numeric; species tolerances, one per taxon.}

  \item{h}{numeric; species abundance at optima, one per
  taxon.}
}
\value{
a matrix of simulated counts.
}
\description{
Simulates species abundances along a single gradient
simulated from a Poisson distribution assuming a Gaussian
response model.
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
Opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
Tol <- rep(0.25, 5)
H <- rep(20, 5)
y <- sim1dPoisson(x1, Opt, Tol, H)
}
\author{
Gavin L. Simpson
}

