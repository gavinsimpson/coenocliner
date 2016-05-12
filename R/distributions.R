##' @title Wrappers to random number generators for use with coenocliner
##'
##' @description These functions are simple wrappers around existing random number generators in R to provide stochastic count data for simulated species.
##'
##' @references Bolker, B.M. (2008) \emph{Ecological Models and Data
##' in R.} Princeton University Press.
##'
##' @param n the number of random draws, equal to number of species times the number of gradient locations.
##' @param mu the mean or expectation of the distribution. For \code{Bernoulli}, \code{Binomial}, and \code{BetaBinomial()} this is the probability of occurrence as given by the response function.
##' @param alpha numeric; dispersion parameter for the negative binomial distribution. May be a vector of length \code{length(mu)}. The NB2 parametrization of the negative binomial is used here, in which \eqn{\alpha} is positively related to the amount of extra dispersion in the simulated data. As such, where \eqn{\alpha = 0}, we would have a Poisson distribution. \code{alpha} can be supplied a value of \code{0}, in which case \code{NegBin} and \code{ZINB} return random draws from the Poisson or zero-inflated Poisson distributions, respectively. Negative values of \code{alpha} are not allowed and will generate an error.
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
    if (any(alpha < 0L)) {
        stop("Negative values of 'alpha' are not supported")
    }
    alpha <- rep(alpha, length.out = length(mu))
    ind <- alpha > 0L
    mu[ind] <- mu[ind] * rgamma(sum(ind), shape = 1/alpha[ind], rate = 1/alpha[ind])
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
##' @importFrom stats rbinom
##'
##' @param size numeric; binomial denominator, the total number of individuals counted for example
`Binomial` <- function(n, mu, size) {
    rbinom(n = n, size = size, prob = mu)
}

##' @rdname distributions
##'
##' @param theta numeric; a positive \emph{inverse} overdispersion parameter for the Beta-Binomial distribution. Low values give high overdispersion. The variance is  \code{size*mu*(1-mu)*(1+(size-1)/(theta+1))} (Bolker, 2008)
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
##' @param zprobs numeric; zero-inflation parameter giving the proportion of extraneous zeros. Must be in range \eqn{0 \dots 1}{0 to 1}.
##' @importFrom stats runif rpois
`ZIP` <- function(n, mu, zprobs) {
    ifelse(runif(n) > zprobs, rpois(n, lambda = mu), 0L)
}

##' @rdname distributions
##'
##' @importFrom stats rpois rgamma runif
`ZINB` <- function(n, mu, alpha, zprobs) {
    if (any(alpha < 0L)) {
        stop("Negative values of 'alpha' are not supported")
    }
    alpha <- rep(alpha, length.out = length(mu))
    ind <- alpha > 0L
    mu[ind] <- mu[ind] * rgamma(sum(ind), shape = 1/alpha[ind], rate = 1/alpha[ind])
    ifelse(runif(n) > zprobs,
           rpois(n, lambda = mu),
           0L)
}

##' @rdname distributions
##'
##' @importFrom stats rbinom runif
## Zero-inflated Binomial
`ZIB` <- function(n, mu, size, zprobs) {
    ifelse(runif(n) > zprobs,
           rbinom(n, size = size, prob = mu),
           0L)
}

##' @rdname distributions
##'
##' @importFrom stats runif
## Zero-inflated Beta-Binomial
`ZIBB` <- function(n, mu, size, theta, zprobs) {
    ifelse(runif(n) > zprobs,
           BetaBinomial(n, mu = mu, size = size, theta = theta),
           0L)
}
