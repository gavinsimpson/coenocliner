`NegBin` <- function(n, mu, alpha) {
    mu <- mu * rgamma(n, shape = alpha, rate = 1/alpha)
    rpois(n, lambda = mu)
}

`Poisson` <- function(n, mu) {
    rpois(n, lambda = mu)
}

`Binomial` <- function(n, mu) {
    rbinom(n, size = 1, prob = mu)
}
