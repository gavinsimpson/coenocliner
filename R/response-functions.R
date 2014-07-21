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

`Beta` <- function(x1, x2 = NULL,
                   m1, A01, r1, alpha1, gamma1,
                   m2 = NULL, A02 = NULL, r2 = NULL,
                   alpha2 = NULL, gamma2 = NULL) {
    sim <- if (is.null(x2)) {
        ## check alpha and gamma are positive as this gives unimodal curves
        stopifnot(alpha1 > 0)
        stopifnot(gamma1 > 0)

        ## number of gradient locations
        ## n <- length(x1)

        ## repeat parameters to length of x1
        ## alpha1 <- rep(alpha1, length.out = n)
        ## gamma1 <- rep(gamma1, length.out = n)
        ## m1 <- rep(m1, length.out = n)
        ## r1 <- rep(r1, length.out = n)
        ## A01 <- rep(A01, length.out = n)

        ## store parameters
        params <- list(x1 = x1, m1 = m1, A01 = A01, r1 = r1,
                       alpha1 = alpha1,
                       gamma1 = gamma1)

        ## some derived parameters, eqns 3 and 4 from Minchin 1987
        b <- alpha1 / (alpha1 + gamma1)
        d <- b^alpha1 * (1 - b)^gamma1

        ## ranges of spp on gradient
        lwr <- m1 - (r1 * b)
        upr <- m1 + (r1 * (1 - b))

        ## output vector A of abundances
        ## set to zero as outside ranges above, this is the abundance
        ## A <- numeric(length = n)
        A <- numeric(length = length(x1))

        ## observations where grad locations x within range of spp
        ind <- x1 >= lwr & x1 <= upr

        ## update parameters - simplifies function call for beta response
        x1 <- x1[ind]
        alpha1 <- alpha1[ind]
        gamma1 <- gamma1[ind]
        m1 <- m1[ind]
        r1 <- r1[ind]
        A01 <- A01[ind]
        b <- b[ind]
        d <- d[ind]

        ## update A with beta response values
        xmr <- (x1 - m1) / r1
        term1 <- xmr + b
        term2 <- 1 - (xmr + b)
        A[ind] <- (A01 / d) * (term1)^alpha1 *
            (term2)^gamma1 ## eqn (2) in Minchin

        ## return
        attr(A, "params") <- params
        A
    } else {
        stop("Bivariate generalised beta response not yet implemented")
    }
    sim
}
