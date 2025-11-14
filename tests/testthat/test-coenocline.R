## Tests for main coenocline() function

## Load packages
library("testthat")
library("coenocliner")

## set up for tests
x <- seq(from = 4, to = 6, length = 100)
opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
tol <- rep(0.25, 5)
h <- rep(20, 5)

## simulate
set.seed(1)
sim <- coenocline(
  x,
  responseModel = "gaussian",
  params = cbind(opt = opt, tol = tol, h = h),
  countModel = "poisson"
)

test_that("coenocline() returns an integer matrix", {
  expect_s3_class(sim, "coenocline")
  expect_s3_class(sim, "matrix")
  expect_type(sim, "integer")
})

test_that("coenocline() returns matrix with correct number of species", {
  expect_identical(NCOL(sim), length(opt))
})

test_that("coenocline() returns matrix with correct number of samples (rows)", {
  expect_identical(NROW(sim), length(x))
})

## simulate
x <- seq(from = 4, to = 6, length = 100)
y <- seq(from = 1, to = 100, length = 100)
optx <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
opty <- c(5, 50, 75, 10, 60)
tolx <- rep(0.25, 5)
toly <- rep(2, 5)
h <- rep(30, 5)

set.seed(1)
sim <- coenocline(cbind(x, y),
                  responseModel = "gaussian",
                  params = list(px = cbind(opt = optx, tol = tolx, h = h),
                                py = cbind(opt = opty, tol = toly)),
                  countModel = "poisson")

test_that("coenocline() returns an integer matrix with x and y gradients", {
  expect_s3_class(sim, "coenocline")
  expect_s3_class(sim, "matrix")
  expect_type(sim, "integer")
})

test_that("coenocline() returns matrix with correct number of species with x and y gradients", {
  expect_identical(NCOL(sim), length(opt))
})

test_that("coenocline() returns matrix with correct number of samples with x and y gradients", {
  expect_identical(NROW(sim), length(x))
})

test_that("coenocline() works with parameters as lists", {
  lp <- list(px = list(opt = optx, tol = tolx, h = h),
             py = list(opt = opty, tol = toly))
  set.seed(1)
  sim2 <- coenocline(cbind(x, y),
                     responseModel = "gaussian",
                     params = lp,
                     countModel = "poisson")
  expect_s3_class(sim2, "coenocline")
  expect_s3_class(sim2, "matrix")
  expect_identical(NCOL(sim2), length(opt))
  expect_identical(NROW(sim2), length(x))
  expect_identical(sim2, sim)
})

test_that("coenocline() works with gradient values supplied as a list", {
  lp <- list(px = list(opt = optx, tol = tolx, h = h),
             py = list(opt = opty, tol = toly))
  set.seed(1)
  sim3 <- coenocline(list(x = x, y = y),
                     responseModel = "gaussian",
                     params = lp,
                     countModel = "poisson")
  expect_s3_class(sim3, "coenocline")
  expect_s3_class(sim3, "matrix")
  expect_identical(NCOL(sim3), length(opt))
  expect_identical(NROW(sim3), length(x))
  expect_identical(sim3, sim)
})
