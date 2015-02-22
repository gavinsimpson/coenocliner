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
##'
##' @examples
##' ## Poisson counts along a single gradient, Gaussian response
##' ## =========================================================
##'
##' x <- seq(from = 4, to = 6, length = 100)
##' opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
##' tol <- rep(0.25, 5)
##' h <- rep(20, 5)
##'
##' ## simulate
##' set.seed(1)
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson")
##' head(y)
##'
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson",
##'                 expectation = TRUE)
##' plot(y, type = "l", lty = "solid")
`plot.coenocline` <- function(x, type = "p", ...) {
    locs <- locations(x)
    ord <- order(locs)
    matplot(locs[ord], x[ord, ], type = type, ...)
}

##' @export
##' @rdname plot.coenocline
##'
##' @importFrom graphics lines matlines
`lines.coenocline` <- function(x, ...) {
    locs <- locations(x)
    ord <- order(locs)
    matlines(locs[ord], x[ord, ], ...)
}
