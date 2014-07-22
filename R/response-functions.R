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
