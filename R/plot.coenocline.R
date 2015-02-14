##' @title Plot species simulations along gradients
##'
##' @description A simple S3 \code{\link{plot}} method for coenocline simulations.
##'
##' @param x an object of class \code{"coenocline"}, the result of a call to \code{\link{coenocline}}.
##' @param type character; the type of plot to produce. See \code{\link{plot.default}} for details.
##' @param ... additional arguments to \code{\link{matplot}}.
##' @return A plot is drawn on the current device.
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @rdname plot.coenocline
##'
##' @keywords hplot
##'
##' @importFrom graphics plot matplot
`plot.coenocline` <- function(x, type = "p", ...) {
    matplot(attr(x, "locations"), x, type = type, ...)
}

##' @export
##' @rdname plot.coenocline
##'
##' @importFrom graphics lines matlines
`lines.coenocline` <- function(x, ...) {
    matlines(attr(x, "locations"), x, ...)
}
