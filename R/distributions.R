##' @title Wrappers to random number generators for use with coenocliner
##'
##' @param n the number of random draws, equal to number of species times the number of gradient locations.
##' @param mu the mean or expectation of the distribution. For \code{Binomial} this is the probability of occurrence as given by the response function.
##' @param alpha numeric; parameter for the negative binomial distribution.
##'
##' @return a vector of random draws from the stated distribution.
##'
##' @author Gavin L. Simpson
##'
##' @rdname distributions
##'
##' @name distributions
`NegBin` <- function(n, mu, alpha) {
    mu <- mu * rgamma(n, shape = alpha, rate = 1/alpha)
    rpois(n, lambda = mu)
}

##' @rdname distributions
`Poisson` <- function(n, mu) {
    rpois(n, lambda = mu)
}

##' @rdname distributions
`Binomial` <- function(n, mu) {
    rbinom(n, size = 1, prob = mu)
}
