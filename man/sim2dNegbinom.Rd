\name{sim2dNegbinom}
\alias{sim2dNegbinom}
\title{Simulate negative binomial counts along two, possibly
correlated gradients}
\usage{
sim2dNegbinom(x1, x2, opt1, tol1, h, opt2, tol2, corr = 0, alpha)
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

  \item{alpha}{numeric; dispersion parameter for the
  negative binomial.}
}
\value{
A matrix of simulate counts.
}
\description{
Simulates species abundances along two possibly correlated
gradients via simulation from a negative binomial
distribution assuming a Gaussian response model.
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
Opt1 <- seq(4, 6, length = 5)
Tol1 <- rep(0.25, 5)
x2 <- runif(300, min = 2, max = 20)
Opt2 <- seq(2, 20, length = 5)
Tol2 <- rep(1, 5)
H <- rep(20, 5)
y <- sim2dNegbinom(x1, x2, Opt1, Opt2, Tol1, Tol2, H,
                   corr = 0.5, alpha = 1.1)
}
\author{
Gavin L. Simpson
}

