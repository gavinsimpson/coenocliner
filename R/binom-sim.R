##' Simulates species occurrence along a single gradient 
##' from a binomial distribution assuming a Gaussian response model.
##'
##' The response model used is the classical Gaussian response with
##' parameters \eqn{\mu}{mu}, the species optimum, \eqn{t}, the
##' species tolerance, and \eqn{h}, the height of the response curve
##' at the optimum (see \code{\link{gaussianResponse}}).
##' 
##' If \code{expectation = FALSE} species occurrence (either presence or absence) 
##' is drawn randomly from a binomial distribution with probability given
##' by the Gaussian response.
##' If \code{expectation = TRUE} then a continuous probability of presence (according to 
##' the Gaussian response) is returned.
##'
##' @section Note:
##' When called with \code{expectation = FALSE} the function does not use
##' the pseudorandom number generator, but when called with the defaults
##' a single call to pseudorandom number generator is made.
##'
##' @title Simulate species occurrence along a single gradient
##' @param x numeric; gradient locations.
##' @param opt numeric; species optima, one per taxon.
##' @param tol numeric; species tolerances, one per taxon.
##' @param h numeric between 0 and 1; species probability of presence at optima, one per taxon.
##' @param expectation logical; should expectations (probability of presence) rather than
##' discrete occurrence be returned?
##'
##' @return a matrix of simulated occurrences (if \code{expectation = FALSE}).
##' A matrix of probability of presence along the environmental gradient otherwise.
##'
##' @author F. Rodriguez-Sanchez, based on code by Gavin L. Simpson
##'
##' @export
##' @seealso \code{\link{gaussianResponse}}
##'
##' @importFrom stats rbinom
##'
##' @examples
##' ## One species:
##' set.seed(1)
##' x <- seq(0, 30, length.out = 100)
##' Opt <- 18
##' Tol <- 3
##' H <- 0.7
##' y <- sim1dBinom(x, Opt, Tol, H, expectation = TRUE)
##' plot(x, y, main = "Probability of Presence", type="l")
##' 
##' y <- sim1dBinom(x, Opt, Tol, H, expectation = FALSE)
##' plot(x, y, main = "Occurrence")
##' 
##' 
##' ## Multiple species
##' set.seed(1)
##' x1 <- seq(from = 4, to = 6, length = 100)
##' Opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
##' Tol <- rep(0.25, 5)
##' H <- rep(0.7, 5)
##' y <- sim1dBinom(x1, Opt, Tol, H)
##'
##' y <- sim1dBinom(x1, Opt, Tol, H, expectation = TRUE)
##' matplot(x1, y, type = "l", lty = "solid")
##'
`sim1dBinom` <- function(x, opt, tol, h, expectation = FALSE) {
    n <- length(x)
    if (any(h < 0 | h > 1)) stop("h must be a probability between 0 and 1")
    ex <- expandGauss(x, opt, tol, h)
    mu <- gaussianResponse(x = ex[,"x"], opt = ex[,"opt"],
                           tol = ex[,"tol"], h = ex[,"h"])
    if (expectation) {
        sim <- mu
    } else {
        sim <- rbinom(nrow(ex), 1, mu)
    }
    sim <- matrix(sim, nrow = n)
    sim
}

##' Simulates species occurrence along two possibly correlated gradients
##' via simulation from a binomial distribution assuming a
##' Gaussian response model.
##'
##' The response model used is the classical Gaussian response with
##' parameters \eqn{\mu}{mu}, the species optimum, \eqn{t}, the
##' species tolerance, and \eqn{h}, the height of the response curve
##' at the optimum (see \code{\link{biGaussianResponse}}).
##'
##' If \code{expectation = FALSE} species occurrence (either presence or absence) 
##' is drawn randomly from a binomial distribution with probability given
##' by the Gaussian response.
##' If \code{expectation = TRUE} then a continuous probability of presence (according to 
##' the Gaussian response) is returned.
##'
##' @section Note:
##' When called with \code{expectation = FALSE} the function does not use
##' the pseudorandom number generator, but when called with the defaults
##' a single call to pseudorandom number generator is made.
##'
##' @title Simulate species occurrence along two, possibly correlated gradients
##' @param x1 numeric; locations along gradient 1.
##' @param x2 numeric; locations along gradient 2.
##' @param opt1 numeric; species optima on gradient 1, one per taxon.
##' @param tol1 numeric; species tolerance on gradient 1, one per taxon.
##' @param h numeric between 0 and 1; species probability of presence at the optimum.
##' @param opt2 numeric; species optima on gradient 2, one per taxon.
##' @param tol2 numeric; species tolerance on gradient 2, one per taxon.
##' @param corr numeric; the correlation between gradients.
##' @param expectation logical; should expectations (probability of presence) rather than
##' discrete occurrence be returned?
##'
##' @return a matrix of simulated occurrences (if \code{expectation = FALSE}).
##' A matrix of probability of presence along the environmental gradients otherwise.
##'
##' @author F. Rodriguez-Sanchez, based on code by Gavin L. Simpson
##'
##' @export
##' @seealso \code{\link{biGaussianResponse}}
##'
##' @importFrom stats rbinom
##'
##' @examples
##' ## One species:
##' set.seed(1)
##' x1 <- runif(100, 0, 30)
##' Opt1 <- 18
##' Tol1 <- 5
##' x2 <- runif(100, 400, 1000)
##' Opt2 <- 800
##' Tol2 <- 100
##' H <- 0.8
##' y <- sim2dBinom(x1, x2, Opt1, Tol1, H, Opt2, Tol2, corr = 0, expectation = TRUE)
##'   
##' ## Multiple species
##' set.seed(1)
##' nsp <- 5
##' x1 <- runif(100, min = 0, max = 30)
##' Opt1 <- seq(5, 25, length = nsp)
##' Tol1 <- rep(5, nsp)
##' x2 <- runif(100, min = 400, max = 1000)
##' Opt2 <- seq(500, 900, length = nsp)
##' Tol2 <- rep(100, nsp)
##' H <- rep(0.9, 5)
##' y <- sim2dBinom(x1, x2, Opt1, Tol1, H, Opt2, Tol2, corr = 0.2)
##'
`sim2dBinom` <- function(x1, x2, opt1, tol1, h, opt2, tol2,
                           corr = 0, expectation = FALSE) {
    stopifnot(isTRUE(all.equal(n1 <- length(x1), n2 <- length(x2))))
    if (any(h < 0 | h > 1)) stop("h must be a probability between 0 and 1")
    ex1 <- expandGauss(x1, opt1, tol1, h)
    ex2 <- expandGauss(x2, opt2, tol2, h)[, -4] ## don't need extra h
    colnames(ex1) <- c("x1", "opt1", "tol1", "h")
    colnames(ex2) <- c("x2", "opt2", "tol2")

    mu <- biGaussianResponse(x1 = ex1[, "x1"], x2 = ex2[, "x2"],
                             opt1 = ex1[, "opt1"], tol1 = ex1[, "tol1"],
                             opt2 = ex2[, "opt2"], tol2 = ex2[, "tol2"],
                             h = ex1[, "h"], corr = corr)

    if (expectation) {
        sim <- mu
    } else {
        sim <- rbinom(n1*length(opt1), 1, mu)
    }
    sim <- matrix(sim, nrow = n1)
    sim
}
