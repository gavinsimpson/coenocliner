##' @title Extract the Matrix of Simulated Counts
##'
##' @description Extracts the matrix of simulated counts from a suitable object.
##' @param x an object of class \code{"coenocline"}, the result of a call to \code{\link{coenocline}}.
##' @param ... additional arguments passed to other methods. Currently ignored as the function is not and S3 generic.
##' @return A matrix of simulated counts
##' @author Gavin L. Simpson
##' @export
`simulations` <- function(x, ...) {
    x$simulations
}
