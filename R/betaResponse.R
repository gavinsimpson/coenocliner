##' @title Generalized beta function of Minchin (1987)
##'
##' @description Species abundances are generated from a generalized beta
##' function. The simple beta function is generalized to any range of
##' \code{x}, the gradient locations. Parameters specify the position
##' of the modal response (\code{m}), the maximum abundance at the mode
##' (\code{A0}), and the range of occurence of the species on the
##' gradient (\code{r}). In addition, two positive shape parameters,
##' \code{alpha} and \code{gamma} constrain the shape of the response.
##'
##' @param x numeric; gradient locations.
##' @param m numeric; location on \code{x} of the peak (modal) abundance.
##' @param A0 numeric; maximum abundance at the mode \code{m}.
##' @param r numeric; range of occurence of species on gradient \code{x}.
##' @param alpha numeric; shape parameter of generalised bete function.
##' Must be positive.
##' @param gamma numeric; shape parameter of generalised bete function.
##' Must be positive.
##'
##' @return A vector of abundance values for the stated parameters. The
##' vector has attribute \code{params}, which contains the specified
##' parameter values.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @references Minchin P.R. (1987) Simulation of multidimensional
##' community patterns: towards a comprehensive model. \emph{Vegetatio}
##' \strong{71}, 145--156.
##'
##' @examples
##' N <- 100
##' X <- seq(from = 3, to = 7, length.out = N)
##' sim <- betaResponse(X, m = 5, A0 = 50, r = 2.5, alpha = 1.5,
##'                     gamma = 0.5)
##' plot(sim ~ X, type = "l")
`betaResponse` <- function(x, m, A0, r, alpha, gamma) {
    ## check alpha and gamma are positive as this gives unimodal curves
    stopifnot(alpha > 0)
    stopifnot(gamma > 0)

    ## number of gradient locations
    n <- length(x)

    ## repeat parameters to length of x
    alpha <- rep(alpha, length.out = n)
    gamma <- rep(gamma, length.out = n)
    m <- rep(m, length.out = n)
    r <- rep(r, length.out = n)
    A0 <- rep(A0, length.out = n)

    ## store parameters
    params <- list(x = x, m = m, A0 = A0, r = r, alpha = alpha, gamma = gamma)

    ## some derived parameters, eqns 3 and 4 from Minchin 1987
    b <- alpha / (alpha + gamma)
    d <- b^alpha * (1 - b)^gamma

    ## ranges of spp on gradient
    lwr <- m - (r * b)
    upr <- m + (r * (1 - b))

    ## output vector A of abundances
    ## set to zero as outside ranges above, this is the abundance
    A <- numeric(length = n)

    ## observations where grad locations x within range of spp
    ind <- x >= lwr & x <= upr

    ## update parameters - simplifies function call for beta response
    x <- x[ind]
    alpha <- alpha[ind]
    gamma <- gamma[ind]
    m <- m[ind]
    r <- r[ind]
    A0 <- A0[ind]
    b <- b[ind]
    d <- d[ind]

    ## update A with beta response values
    xmr <- (x - m) / r
    term1 <- xmr + b
    term2 <- 1 - (xmr + b)
    A[ind] <- (A0 / d) * (term1)^alpha * (term2)^gamma ## eqn (2) in Minchin

    ## return
    attr(A, "params") <- params
    A
}
