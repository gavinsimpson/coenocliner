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
##' @return a matrix of Generalised Beta repsonse model parameters for
##' each gradient location. Specifically, a matrix with columns \code{x},
##' \code{m}, \code{A0}, \code{r}, \code{alpha}, and \code{gamma}.
##'
##' @author Gavin L. Simpson
##'
##' @export
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
