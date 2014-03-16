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
