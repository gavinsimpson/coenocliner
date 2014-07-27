##' @title Wrappers to random number generators for use with coenocliner
##'
##' @param n the number of random draws, equal to number of species times the number of gradient locations.
##' @param mu the mean or expectation of the distribution. For \code{Bernoulli}, \code{Binomial}, and \code{BetaBinomial()} this is the probability of occurrence as given by the response function.
##' @param alpha numeric; parameter for the negative binomial distribution.
##'
##' @return a vector of random draws from the stated distribution.
##'
##' @author Gavin L. Simpson
##'
##' @rdname distributions
##'
##' @name distributions
##'
##' @keywords distribution
##'
##' @importFrom stats rpois
`NegBin` <- function(n, mu, alpha) {
    mu <- mu * rgamma(n, shape = alpha, rate = 1/alpha)
    rpois(n, lambda = mu)
}

##' @rdname distributions
`Poisson` <- function(n, mu) {
    rpois(n, lambda = mu)
}

##' @rdname distributions
##'
##' @importFrom stats rbinom
`Bernoulli` <- function(n, mu) {
    rbinom(n, size = 1, prob = mu)
}

##' @rdname distributions
##'
##' @param size numeric; binomial denominator, the total number of individuals counted for example
`Binomial` <- function(n, mu, size) {
    rbinom(n = n, size = size, prob = mu)
}

##' @rdname distributions
##'
##' @param tau numeric; the overdispersion parameter for the Beta-Binomial distribution. This is actually tau^2
##'
##' @importFrom stats rbeta
`BetaBinomial` <-  function(n, mu, size, tau) {
    ## tau is tau^2
    t1 <- tau * size
    t2 <- t1 - size - tau + 1
    t3 <- 1 / (1 + t1 - tau)
    t4 <- t2 * t3
    a  <- -t4 * mu
    b  <- t4 * (mu - 1)
    mu <- rbeta(n = n, shape1 = a, shape2 = b)
    rbinom(n = n, size = size, prob = mu)
}
