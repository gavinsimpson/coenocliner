`Gaussian` <- function(x1, x2 = NULL,
                       opt1, tol1, opt2 = NULL,
                       tol2 = NULL, h, corr = NULL) {
    sim <- if (is.null(x2)) {
        h * exp(-((x1 - opt1)^2/(2 * tol1^2)))
    } else {
        if (is.null(corr)) {
            corr <- 0
        }
        h * exp(-(1/(2 * (1 - corr^2))) * (((x1 - opt1)/tol1)^2 +
            ((x2 - opt2)/tol2)^2 - ((2 * corr) * ((x1 - opt1)/tol1) *
            ((x2 - opt2)/tol2))))
    }
    sim
}
