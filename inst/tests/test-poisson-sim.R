## Tests for Possion count simulation functions

## Load packages
library("testthat")
library_if_available("coenocliner")

context("Testing Possion count simulation")

### 1 gradient

## Set up options for simulations
set.seed(1)
x1 <- seq(from = 4, to = 6, length = 100)
Opt <- seq(4, 6, length = 5)
Tol <- rep(0.25, 5)
H <- rep(20, 5)
y <- sim1dPoisson(x1, Opt, Tol, H)

test_that("sim1dPoisson returns correctly a matrix", {
    expect_that(y, is_a("matrix"))
    expect_that(NROW(y), equals(length(x1)))
    expect_that(NCOL(y), equals(length(Opt)))
})

### 2 correlated gradients

## Set up options for simulations
set.seed(1)
x1 <- runif(300, min = 4, max = 6)
Opt1 <- seq(4, 6, length = 5)
Tol1 <- rep(0.25, 5)
x2 <- runif(300, min = 2, max = 20)
Opt2 <- seq(2, 20, length = 5)
Tol2 <- rep(1, 5)
H <- rep(20, 5)
y <- sim2dPoisson(x1, x2, Opt1, Opt2, Tol1, Tol2, H,
                  corr = 0.5)

test_that("sim2dPoisson returns correctly a matrix", {
    expect_that(y, is_a("matrix"))
    expect_that(NROW(y), equals(length(x1)))
    expect_that(NCOL(y), equals(length(Opt1)))
})
