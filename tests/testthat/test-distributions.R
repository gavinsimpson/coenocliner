## Tests for response model functions

## Load packages
library("testthat")
library("coenocliner")

context("Testing RNG wrappers")

n <- 10L

test_that("NegBin returns the right number of random draws", {
    expect_identical(length(NegBin(n, mu = 2, alpha = 1)), n)
})

test_that("Poisson returns the right number of random draws", {
    expect_identical(length(Poisson(n, mu = 2)), n)
})

test_that("Bernoulli returns the right number of random draws", {
    expect_identical(length(Bernoulli(n, mu = 0.5)), n)
})

test_that("Binomial returns the right number of random draws", {
    expect_identical(length(Binomial(n, mu = 0.5, size = 10)), n)
})

test_that("BetaBinomial returns the right number of random draws", {
    expect_identical(length(BetaBinomial(n, mu = 0.5, size = 10,
                                         theta = 2)), n)
})

test_that("ZIP returns the right number of random draws", {
    expect_identical(length(ZIP(n, mu = 2, zprobs = 0.5)), n)
})

test_that("ZINB returns the right number of random draws", {
    expect_identical(length(ZINB(n, mu = 2, alpha = 1, zprobs = 0.5)), n)
})

test_that("ZIB returns the right number of random draws", {
    expect_identical(length(ZIB(n, mu = 0.5, size = 10, zprobs = 0.5)), n)
})

test_that("ZIBB returns the right number of random draws", {
    expect_identical(length(ZIBB(n, mu = 0.5, size = 10,
                                 theta = 2, zprobs = 0.5)), n)
})

test_that("NegBin handles negative alpha correctly", {
    expect_error(NegBin(n, mu = 2, alpha = -1),
                 regexp = "Negative values of 'alpha' are not supported")
})

test_that("NegBin returns same as Poisson for alpha = 0", {
    mu <- 2
    set.seed(1)
    rnd1 <- NegBin(n, mu = mu, alpha = 0)
    set.seed(1)
    rnd2 <- Poisson(n, mu = mu)
    expect_identical(rnd1, rnd2)
})

test_that("ZINB handles negative alpha correctly", {
    expect_error(ZINB(n, mu = 2, alpha = -1, zprobs = 0.2),
                 regexp = "Negative values of 'alpha' are not supported")
})

test_that("ZINB returns same as ZIP for alpha = 0", {
    mu <- 2
    zprobs <- 0.2
    set.seed(1)
    rnd1 <- ZINB(n, mu = mu, alpha = 0, zprobs = zprobs)
    set.seed(1)
    rnd2 <- ZIP(n, mu = mu, zprobs = zprobs)
    expect_identical(rnd1, rnd2)
})

test_that("NegBin and ZINB work with vector alpha argument", {
    mu <- 2
    zprobs <- 0.2
    set.seed(1)
    rndNegBin <- NegBin(n, mu = mu, alpha = 0)
    rndZINB   <- ZINB(n, mu = mu, alpha = 0, zprobs = zprobs)
    expect_is(rndNegBin, "integer")
    expect_length(rndNegBin, 10L)
    expect_is(rndZINB, "integer")
    expect_length(rndZINB, 10L)
})
