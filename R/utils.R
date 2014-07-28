##' @title List parameters of species response models
##'
##' @description Returns the parameters of the indicated response model.
##'
##' @param model character; the species response model for which parameters will be listed
##'
##' @return A character vector of parameters. The species response model is returned as attribute \code{"model"}. Attribute \code{"onlyx"} is a logical vector indicating which, if any, of the parameters are intended to be supplied only once per species and not for both gradients.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @seealso \code{\link{Gaussian}} and \code{\link{Beta}} for the species response model functions themselves.
##'
##' @examples
##' showParams("gaussian")
`showParams` <- function(model = c("gaussian", "beta")) {
    if (missing(model)) {
        stop("The species response model be provided.")
    }
    model <- match.arg(model)
    switch(model,
           gaussian = {
               params <- c("opt", "tol", "h")
               attr(params, "onlyx") <- c(rep(FALSE, 2), TRUE)
               attr(params, "model") <- "Gaussian"
           },
           beta = {
               params <- c("m", "r", "alpha", "gamma", "A0")
               attr(params, "onlyx") <- c(rep(FALSE, 4),TRUE)
               attr(params, "model") <- "Generalise Beta"
           })
    class(params) <- "modelParams"
    params
}

##' @export
`print.modelParams` <- function(x, ...) {
    onlyx <- attr(x, "onlyx")
    if (length(x)) {
        msg <- paste("Species response model:", attr(x, "model"))
        writeLines(strwrap(msg))
        writeLines(strwrap("Parameters:", prefix = "\n"))
        print(paste0(x, ifelse(onlyx, "*", "")), quote = FALSE, ...)
        if (any(onlyx)) {
            msg <- "Parameters marked with '*' are only supplied once"
            writeLines(strwrap(msg, prefix = "\n\t"))
        }
    }
    invisible(x)
}
