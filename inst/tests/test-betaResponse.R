## Tests for Beta response function

## Load packages
library("testthat")
library_if_available("coenocliner")

context("Testing Beta response function")


## Set up options for simulations
N <- 100
X <- seq(from = 3, to = 7, length.out = N)
sim <- betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = 1.5,
                    gamma = 0.5)

test_that("betaResponse returns a vector", {
    expect_that(sim, is_a("numeric"))
})

test_that("betaResponse fails with non-positive shape parameters", {
    expect_that(betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = -1.5,
                             gamma = 0.5),
                throws_error())
    expect_that(betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = 1.5,
                             gamma = -0.5),
                throws_error())
    expect_that(betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = -1.5,
                             gamma = -0.5),
                throws_error())
})
