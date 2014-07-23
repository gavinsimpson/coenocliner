##' @title Species response models for coenocline simulation
##'
##' @description Parameterise species response curves along one or two gradients according to a Gaussian or generalised beta response model.
##'
##' @details \code{Gaussian()} and \code{Beta()} return values from appropriately parameterised Gaussian or generalised beta response models respectively. Parameters for the primary (\code{x}) and secondary (\code{y}) gradients are supplied as lists via arguments \code{px} and \code{py}. Parameters are supplied in the form of vectors, one per parameter. These vectors must be supplied to named components in the respective lists. The names of the components must match the parameters of the required response model.
##'
##' For \code{Gaussian()} the following named components must be supplied:
##' \describe{
##'   \item{opt}{the species optima}
##'   \item{tol}{the species tolerances}
##'   \item{h}{the heights of the response curves at the optima. This parameter should only be supplied to \code{px}; in the case of simulations along two gradients, the height of the response curve applies to both gradients and is the hieght of a bivariate Guassian distribution at the bivariate optima.}
##' }
##'
##' For \code{Beta()} the following named components must be supplied:
##' \describe{
##'   \item{A0}{The heights of the species response curves at their modes. Like the parameter \code{h} for the Gaussian response, this parameter should only be passed via \code{px}; in the case of simulations along two gradients, the height of the response curve applies to both gradients and is the height of a bivariate generalised beta distribution at the bivariate mode.}
##'   \item{m}{the locations on the gradient of the modal abundance (the species optima)}
##'   \item{r}{the ranges of occurrence of species on the gradient}
##'   \item{alpha}{a shape parameter. With \code{gamma}, \code{alpha} informs the shape of the response curve and control the skewness and kurtosis of the curve. Only positive values are allowed, which lead to unimodal response curves. If \code{alpha} is equal to \code{gamma}, the species response curve is symmetric, otherwise an asymmetric curve is generated.}
##'   \item{gamma}{a shape parameter. With \code{alpha}, \code{gamma} informs the shape of the response curve and control the skewness and kurtosis of the curve. Only positive values are allowed, which lead to unimodal response curves. If \code{gamma} is equal to \code{alpha}, the species response curve is symmetric, otherwise an asymmetric curve is generated.}
##' }
##'
##' See the examples here and in \code{\link{coenocline}} for details on how to set up calls to these species response functions.
##'
##' @param x numeric; locations of observations on the primary gradient.
##' @param y numeric; locations of observations on the secondary gradient. Can be missing is only a single gradient is required.
##' @param px a list of named elements, each of which is a vector of numeric parameter values for the species response on the primary gradient \code{x}. See Details for further information on the required parameters.
##' @param py a list of named elements, each of which is a vector of numeric parameter values for the species response on the secondary gradient \code{y}. See Details for further information on the required parameters.
##' @param corr numeric; the correlation between gradients \code{x} and \code{y}. Only applies to \code{Gaussian()}.
##'
##' @return A numeric vector of species "abundances" of length equal to \code{length(x)}.
##'
##' @section Note:
##' Currently, the generalised beta distribution is curently implemented for a single gradient.
##'
##' @author Gavin L. Simpson
##'
##' @rdname species-response
##' @name species-response
##' @export
`Gaussian` <- function(x, y = NULL, px, py = NULL, corr = NULL) {
    pars <- c("opt", "tol", "h")
    sim <- if (is.null(y)) {
        stopifnot(length(px) == 3L)
        check <- pars %in% names(px)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'px'. Check names."))
        }
        if (length(unique(sapply(px, length))) != 1L) {
            stop("Parameter vectors supplied in 'px' are of differing lengths.")
        }

        ## Compute Gaussian response
        px[["h"]] * exp(-((x - px[["opt"]])^2/(2 * px[["tol"]]^2)))
    } else {
        stopifnot(all.equal(length(x), length(y)))
        stopifnot(length(px) == 3L)
        stopifnot(length(py) == 2L)
        if (is.null(corr)) {
            corr <- 0
        }
        check <- pars %in% names(px)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'px'. Check names."))
        }
        check <- pars[1:2] %in% names(py)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'py'. Check names."))
        }
        if (length(unique(sapply(px, length))) != 1L) {
            stop("Parameter vectors supplied in 'px' are of differing lengths.")
        }
        if (length(unique(sapply(py, length))) != 1L) {
            stop("Parameter vectors supplied in 'py' are of differing lengths.")
        }
        px[["h"]] * exp(-(1/(2 * (1 - corr^2))) *
                        (((x - px[["opt"]])/px[["tol"]])^2 +
                         ((y - py[["opt"]])/py[["tol"]])^2 -
                         ((2 * corr) *
                          ((x - px[["opt"]])/px[["tol"]]) *
                          ((y - py[["opt"]])/py[["tol"]]))))
    }
    sim
}

##' @rdname species-response
##' @export
`Beta` <- function(x, y = NULL, px, py = NULL) {
    pars <- c("A0", "m", "r", "alpha", "gamma")
    sim <- if (is.null(y)) {
        ## checks on parameters
        stopifnot(length(px) == 5L)
        check <- pars %in% names(px)
        if (!all(check)) {
            stop(paste("One or more of", pars, "not in 'px'. Check names."))
        }
        if (length(unique(sapply(px, length))) != 1L) {
            stop("Parameter vectors supplied in 'px' are of differing lengths.")
        }

        ## check alpha and gamma are positive as this gives unimodal curves
        stopifnot(px[["alpha"]] > 0)
        stopifnot(px[["gamma"]] > 0)

        ## store parameters
        ## params <- list(x = x, m = m, A0 = A0, r = r,
        ##                alpha = alpha,
        ##                gamma = gamma)

        ## some derived parameters, eqns 3 and 4 from Minchin 1987
        b <- px[["alpha"]] / (px[["alpha"]] + px[["gamma"]])
        d <- b^px[["alpha"]] * (1 - b)^px[["gamma"]]

        ## ranges of spp on gradient
        lwr <- px[["m"]] - (px[["r"]] * b)
        upr <- px[["m"]] + (px[["r"]] * (1 - b))

        ## output vector A of abundances
        ## set to zero as outside ranges above, this is the abundance
        ## A <- numeric(length = n)
        A <- numeric(length = length(x))

        ## observations where grad locations x within range of spp
        ind <- x >= lwr & x <= upr

        ## update parameters - simplifies function call for beta response
        x <- x[ind]
        px[["alpha"]] <- px[["alpha"]][ind]
        px[["gamma"]] <- px[["gamma"]][ind]
        px[["m"]] <- px[["m"]][ind]
        px[["r"]] <- px[["r"]][ind]
        px[["A0"]] <- px[["A0"]][ind]
        b <- b[ind]
        d <- d[ind]

        ## update A with beta response values
        xmr <- (x - px[["m"]]) / px[["r"]]
        term <- xmr + b
        term2 <- 1 - (xmr + b)
        A[ind] <- (px[["A0"]] / d) * (term)^px[["alpha"]] *
            (term2)^px[["gamma"]] ## eqn (2) in Minchin

        ## return
        attr(A, "params") <- list(px = px, py = py)
        A
    } else {
        stop("Bivariate generalised beta response not yet implemented")
    }
    sim
}
