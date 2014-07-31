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
##' @importFrom stats rpois rgamma
`NegBin` <- function(n, mu, alpha) {
    mu <- mu * rgamma(n, shape = alpha, rate = 1/alpha)
    rpois(n, lambda = mu)
}

##' @rdname distributions
##'
##' @importFrom stats rpois
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
##' @param theta numeric; a positive overdispersion parameter for the Beta-Binomial distribution.
##'
##' @importFrom stats rbeta rbinom
`BetaBinomial` <- function(n, mu, size, theta) {
    ## follows Bolker (2008) and derive a and b from mu and theta
    ## mu == pi (or prob, or p)
    a <- theta * mu
    b <- theta * (1 - mu)
    rbinom(n = n, size = size,
           prob = rbeta(n = n, shape1 = a, shape2 = b))
}

##' @rdname distributions
##'
##' @param gamma numeric; zero-inflation parameter. Leads to the probability of a zero, \eqn{\pi}{pi}, in the binomial part of the ZIP via \eqn{\pi = e^\gamma / (1 + e^\gamma)}{pi = e^gamma / (1 + e^gamma)}. Setting \code{gamma = 0} gives a probability of zero from the binomial part of \eqn{\pi = 0.5}{pi = 0.5}.
##' @importFrom stats rbinom rpois
`ZIP` <- function(n, mu, gamma) {
    e <- exp(1)
    pi <- e^gamma / (1 + e^gamma)
    pres <- rbinom(n, size = 1, prob = pi)
    rand <- ifelse(pres > 0, rpois(n, lambda = mu), 0)
    rand
}

##' @rdname distributions
##'
##' @importFrom stats rpois rgamma rbinom
`ZINB` <- function(n, mu, alpha, gamma) {
    e <- exp(1)
    pi <- e^gamma / (1 + e^gamma)
    pres <- rbinom(n, size = 1, prob = pi)
    rand <- ifelse(pres > 0,
                   rpois(n, lambda = mu * rgamma(n, shape = alpha,
                            rate = 1/alpha)),
                   0)
    rand
}
