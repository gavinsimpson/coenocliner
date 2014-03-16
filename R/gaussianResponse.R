##' @title Univariate Gaussian response model
##'
##' @param x numeric; gradient locations.
##' @param opt numeric; species optima, one per taxon.
##' @param tol numeric; species tolerances, one per taxon.
##' @param h numeric; species abundance at optima, one per taxon.
##'
##' @return A vector of expected responses according to the univariate
##' Gaussian response model.
##'
##' @export
##'
##' @author Gavin Simpson
`gaussianResponse` <- function(x, opt, tol, h) {
    h * exp(-((x - opt)^2 / (2*tol^2)))
}

##' @title Bivariate Gaussian response model
##'
##' @param x1 numeric; locations along gradient 1.
##' @param x2 numeric; locations along gradient 2.
##' @param opt1 numeric; species optima on gradient 1, one per taxon.
##' @param tol1 numeric; species toleranceon gradient 1, one per taxon.
##' @param h numeric; species abundance at the optimum
##' @param opt2 numeric; species optima on gradient 2, one per taxon.
##' @param tol2 numeric; species tolerance on gradient 2, one per taxon.
##' @param corr numeric; the correlation between gradients.
##'
##' @return A vector of expected responses according to the bivariate
##' Gaussian repsonse model.
##'
##' @export
##'
##' @author Gavin Simpson
`biGaussianResponse` <- function(x1, x2, opt1, tol1, opt2, tol2, h,
                                 corr = 0) {
    h * exp(-(1/(2*(1 - corr^2))) * (((x1 - opt1) / tol1)^2 +
                                     ((x2 - opt2) / tol2)^2 -
                                     ((2 * corr) * ((x1 - opt1) / tol1) *
                                      ((x2 - opt2) / tol2))))
}
