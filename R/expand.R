##' An \code{\link{expand.grid}}-like utility function for parameters
##' of the Gaussian response model.
##'
##' The response model used is the classical Gaussian response with
##' parameters \eqn{\mu}{mu}, the species optimum, \eqn{t}, the
##' species tolerance, and \eqn{h}, the height of the response curve
##' at the optimum.
##'
##' \code{expandGauss} performs an \code{\link{expand.grid}}-like
##' operation on the gradient locations and species optima, and then
##' replicates the tolerance and height parameters to match the
##' replication of the optima.
##'
##' @title Expand gradient locations and species parameters of the
##' Gaussian response model
##'
##' @param x numeric; gradient locations.
##' @param opt numeric; species optima, one per taxon.
##' @param tol numeric; species tolerances, one per taxon.
##' @param h numeric; species abundance at optima, one per taxon.
##'
##' @return a matrix of Gaussian repsonse model parameters for each
##' gradient location. Specifically, a matrix with columns \code{x},
##' \code{opt}, \code{tol}, and \code{h}.
##'
##' @author Gavin L. Simpson
##'
##' @export
`expandGauss` <- function(x, opt, tol, h) {
    ## do like expand.grid for x and opt, plus expand tol and h
    nx <- length(x)
    no <- length(opt)
    orep <- nx * no
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep/nx)]
    opt <- opt[rep.int(rep.int(seq_len(no), rep.int(nx, no)), 1L)]
    tol <- tol[rep.int(rep.int(seq_len(no), rep.int(nx, no)), 1L)]
    h   <-   h[rep.int(rep.int(seq_len(no), rep.int(nx, no)), 1L)]
    cbind(x, opt, tol, h)
}

##' @title Expand gradient locations and species parameters of the
##' Generalise Beta response
##'
##' @description An \code{\link{expand.grid}}-like utility function
##' for parameters of the Beta response model.
##'
##' @details The response model used is the generalised Beta response with
##' the parameters: the position of the modal response (\code{m}), the
##' maximum abundance at the mode (\code{A0}), and the range of
##' occurence of the species on the gradient (\code{r}). In addition,
##' two positive shape parameters, \code{alpha} and \code{gamma}
##' constrain the shape of the response.\eqn{\mu}{mu}, the species
##' optimum, \eqn{t}, the species tolerance, and \eqn{h}, the height of
##' the response curve at the optimum.
##'
##' \code{expandBeta} performs an \code{\link{expand.grid}}-like operation
##' on the gradient locations and species modela response, and then
##' replicates the remaining Beta response parameters to match the
##' replication of the modal locations.
##'
##' @param x numeric; gradient locations.
##' @param m numeric; species model abundance, one per taxon.
##' @param A0 numeric; species abundance at \code{m}, one per taxon.
##' @param r numeric; range of \code{x} occupied by each species, one
##' per taxon.
##' @param alpha numeric; shape parameter of the generalised beta function.
##' @param gamma numeric; shape parameter of the generalised beta function.
##'
##' @return a matrix of Generalised Beta response model parameters for
##' each gradient location. Specifically, a matrix with columns \code{x},
##' \code{m}, \code{A0}, \code{r}, \code{alpha}, and \code{gamma}.
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
##' # Recreate Fig. 2 of Minchin (1987)
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' x <- 1:100
##'
##' # expand parameter set
##' pars <- expandBeta(x, m, A0, r, alpha, gamma)
##' head(pars)
##'
##' # Simulate expectation under Generalised Beta response
##' sim <- do.call(betaResponse, data.frame(pars))
##' sim <- matrix(sim, nrow = 100)
##'
##' # plot
##' matplot(sim, type = "l", lty = "solid")
`expandBeta` <- function(x, m, A0, r, alpha, gamma) {
    nx <- length(x)
    nm <- length(m)
    orep <- nx * nm
    x     <-     x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep/nx)]
    m     <-     m[rep.int(rep.int(seq_len(nm), rep.int(nx, nm)), 1L)]
    A0    <-    A0[rep.int(rep.int(seq_len(nm), rep.int(nx, nm)), 1L)]
    r     <-     r[rep.int(rep.int(seq_len(nm), rep.int(nx, nm)), 1L)]
    alpha <- alpha[rep.int(rep.int(seq_len(nm), rep.int(nx, nm)), 1L)]
    gamma <- gamma[rep.int(rep.int(seq_len(nm), rep.int(nx, nm)), 1L)]
    cbind(x, m, A0, r, alpha, gamma)
}

##' @title An \code{expand.grid}-like function that repeats sets of
##' vectors for every value in a reference vector.
##'
##' @description The values of \code{x} are repeated for each combination
##' of elements in the vectors supplied via \code{...}, with the first
##' elements of each vector in \code{...} being taken as a set, the
##' second elements as another set, and so on. \code{x} is repeated for
##' each of these sets.
##'
##' @param x numeric; vector of data points which are to be replicated
##' for each of the sets of vectors supplied to \code{...}.
##' @param ... additional vector arguments to be expanded to the correct
##' length. These are taken to be a set of values to be replicated for
##' each of the elements of \code{x}.
##'
##' @return a matrix of replicated vectors, with column names for \code{x}
##' and named arguments passed as \code{...}.
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
##' # Recreate Fig. 2 of Minchin (1987)
##' # Parameters for each of 6 six species
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' # Gradient locations
##' x <- 1:100
##'
##' # expand parameter set
##' pars <- expand(x, m = m, A0 = A0, r = r, alpha = alpha,
##'                gamma = gamma)
##' head(pars)
`expand` <- function(x, ...) {
    ## if `...` is length 1 & matrix use expandMatrix()
    ## otherwise, expandList
    dots <- list(...)
    out <- if (length(dots) == 1L && is.matrix(dots[[1]])) {
        expandMatrix(x, params = dots[[1]])
    } else {
        args <- vector("list", length = length(dots) + 1)
        args[[1]] <- x
        args[-1] <- dots
        names(args) <- c("x", names(dots))
        do.call("expandList", args)
    }
    out
}

`expandList` <- function(x, ...) {
    expandFun <- function(x, r1, r2) {
        x[rep.int(rep.int(seq_len(r2), rep.int(r1, r2)), 1L)]
    }
    dots <- list(...)
    nams <- names(dots)
    ## nams could be NULL or have zero length values
    namsD <- paste0("Var", seq_along(dots))
    if (is.null(nams)) {
        nams <- namsD
    } else if (any(want <- nzchar(nams))){
        namsD[want] <- nams[want]
    }
    nx <- length(x)
    n1 <- length(dots[[1]])
    orep <- nx * n1
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep / nx)]
    other <- vapply(dots, FUN = expandFun,
                    FUN.VALUE = numeric(length = length(x)),
                    r1 = nx, r2 = n1)
    out <- cbind(x, other)
    colnames(out)[-1] <- namsD
    out
}

`expandMatrix` <- function(x, params) {
    nx <- length(x)
    np <- NROW(params)
    orep <- nx * np
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(1L, nx)), orep / nx)]
    other <- params[rep.int(rep.int(seq_len(np), rep.int(nx, np)), 1L), ]
    out <- cbind(x, other)
    colnames(out) <- c("x", colnames(params))
    out
}
