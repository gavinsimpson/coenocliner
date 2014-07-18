##' Simulates species abundances along a single gradient with counts
##' generated from a generalised Beta response model using Poisson or
##' negative binomial random deviates.
##'
##' @description Species abundances are generated from a generalized beta
##' function. The simple beta function is generalized to any range of
##' \code{x}, the gradient locations. Parameters specify the position
##' of the modal response (\code{m}), the maximum abundance at the mode
##' (\code{A0}), and the range of occurence of the species on the
##' gradient (\code{r}). In addition, two positive shape parameters,
##' \code{alpha} and \code{gamma} constrain the shape of the response.
##'
##' If \code{expectation = TRUE} the mean response for the parameters
##' is returned. If \code{expectation = FALSE} counts are drawn randomly
##' from a Poisson or negative binomial distribution with mean (argument
##' \code{lamba}) given by the generalised Beta response.
##'
##' @section Note:
##' When called with \code{expectation = FALSE} the function does not use
##' the pseudorandom number generator, but when called with the defaults
##' a single call to pseudorandom number generator is made.
##'
##' @title Simulate Poisson counts along a single gradient
##' @param x numeric; gradient locations.
##' @param m numeric; location on \code{x} of the peak (modal) abundance.
##' @param A0 A0 numeric; maximum abundance at the mode \code{m}.
##' @param r numeric; range of occurence of species on gradient \code{x}.
##' @param alpha numeric; shape parameter of generalised beta function.
##' @param gamma numeric; shape parameter of generalised beta function.
##' Must be positive.
##' @param count character; which distribution to simulate counts from?
##' One of \code{"poisson"} or \code{"negbin"}.
##' @param expectation logical; should expectations (mean response) be
##' returned?
##' @param nb.alpha numeric; dispersion parameter for the negative binomial.
##'
##' @return a matrix of simulated counts.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @importFrom stats rpois rgamma
##'
##' @examples
##' ## Recreate Fig. 2 of Minchin (1987)
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' x <- 1:100
##'
##' ## Expectations
##' y <- sim1dBeta(x, m, A0, r, alpha, gamma, expectation = TRUE)
##' matplot(y, type = "l", lty = "solid")
##'
##' ## Poisson, simulated
##' set.seed(1)
##' y <- sim1dBeta(x, m, A0, r, alpha, gamma, count = "poisson")
##' matplot(y, type = "p", lty = "solid")
##'
##' ## Negative Binomial, simulated
##' set.seed(1)
##' y <- sim1dBeta(x, m, A0, r, alpha, gamma, count = "negbin",
##'                nb.alpha = 1)
##' matplot(y, type = "p", lty = "solid")
##'
`sim1dBeta` <- function(x, m, A0, r, alpha, gamma,
                        count = c("poisson", "negbin"),
                        expectation = FALSE,
                        nb.alpha) {
    n <- length(x)
    ex <- expandBeta(x, m, A0, r, alpha, gamma)
    lambda <- do.call(betaResponse, data.frame(ex))
    if (expectation) {
        sim <- lambda
    } else {
        count <- match.arg(count)
        sim <- if (identical(count, "poisson")) {
            rpois(NROW(ex), lambda)
        } else {
            nr <- NROW(ex)
            lambda <- lambda * rgamma(nr, shape = nb.alpha, rate = 1/nb.alpha)
            sim <- rpois(nr, lambda)
        }
    }
    sim <- matrix(sim, nrow = n)
    sim
}
