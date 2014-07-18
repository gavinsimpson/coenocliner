## Test `coenocliner` using the `testthat` package

## Setup
library("testthat")
library_if_available("coenocliner")

## Runs the tests in inst/tests
test_package("coenocliner")
