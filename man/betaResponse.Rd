\name{betaResponse}
\alias{betaResponse}
\title{Generalized beta function of Minchin (1987)}
\usage{
betaResponse(x, m, A0, r, alpha, gamma)
}
\arguments{
  \item{x}{numeric; gradient locations.}

  \item{m}{numeric; location on \code{x} of the peak
  (modal) abundance.}

  \item{A0}{numeric; maximum abundance at the mode
  \code{m}.}

  \item{r}{numeric; range of occurence of species on
  gradient \code{x}.}

  \item{alpha}{numeric; shape parameter of generalised bete
  function. Must be positive.}

  \item{gamma}{numeric; shape parameter of generalised bete
  function. Must be positive.}
}
\value{
A vector of abundance values for the stated parameters. The
vector has attribute \code{params}, which contains the
specified parameter values.
}
\description{
Species abundances are generated from a generalized beta
function. The simple beta function is generalized to any
range of \code{x}, the gradient locations. Parameters
specify the the position of the modal response (\code{m}),
the maximum abundance at the mode (\code{A0}), and the
range of occurence of the species on the gradient
(\code{r}). In addition, two positive shape parameters,
\code{alpha} and \code{gamma} constrain the shape of the
response.
}
\examples{
N <- 100
X <- seq(from = 3, to = 7, length.out = N)
sim <- betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = 1.5,
                    gamma = 0.5)
plot(sim ~ X, type = "l")
}
\author{
Gavin L. Simpson
}
\references{
Minchin P.R. (1987) Simulation of multidimensional
community patterns: towards a comprehensive model.
*Vegetatio* **71**, 145--156.
}

