##' @title Simulate species abundance (counts) or occurrence along one or two
##' gradients
##'
##' @description Simulate species abundance (counts) or occurrence along one or two gradients using well-known ecological response models and random draws from a Poisson, negative binomial or binomial distribution.
##'
##' @details \code{coenocline()} is a generic interface to coenocline simulation allowing for easy extension and a consistent interface to a range of species response models and statistical distributions.
##'
##' Two species response models are currently available; the Gaussian response and the generalized beta response model. Random count or occurrence data can be produced via random draws from a suitable distribution; in which case the values obtained from the specoes response function are used as the expectation of the distribution from which random draws are made.
##'
##' Parameters for each species in the response model are supplied via argument \code{params} and can be provided in one of two ways: i) as a list with named components, each of which is a vector containing values for a single parameter for each species, or ii) as a matrix where each column contains the values for a single parameter and the rows represent species. In each case, the names of the list components or the column names of the matrix must be named for the arguments of the function implementing the species distribution model. See the examples.
##'
##' Some species response models may require additional parameters not specified at the per species level. An example is the correlation between gradients in the bivariate Gaussian response model. Such parameters are passed via list \code{extraParams} and must be named accordingly so that they are passed to the corrct argument in the species response function.
##'
##' The species response model defines the mean of expected response. (In the case of a species occurrence, the probability of occurrence is the expectation.) These represent paramterterised distributions. Random count or occurence data can be produced from these distributions by simulation from those distributions. In this case, a count or Bernoulli model is used and random draws from the distribution are made. Currently, the Poisson, the negative binomial, and the Bernoulli distribution are available. Some distributions may need additional parameters beyond the expectation; an example is the \eqn{\alpha}{alpha} parameter of (one parameterisation of) the negative binomial distribution. These parameters are specied via the list \code{countParams}.
##'
##' @param x one of a numeric vector, a list with two components, each a numeric vector, or a matrix with two columns. The vectors are the locations along the gradient(s) at which species responses are to be simulated.
##' @param responseModel character; which species response model to use.
##' @param params a list of vectors each of which are parameters for the response model for each species. Alternatively, a matrix with one column per parameter and a row for each species.
##' @param extraParams a list containing additional parameters required for the response model. Examples include the correlation between gradients in the bivariate Gaussian response model. Components need to be named.
##' @param countModel character; if \code{expectation} is \code{FALSE}, the default, counts (occurrence) are generated using random deviates from the specified distribution.
##' @param countParams a list of additional parameters required to specify the distribution. An example is the parameter \eqn{\alpha}{alpha} in the negative binomial distribution. Components need to be named.
##' @param expectation logical; should the expectation (mean) response be returned (\code{TRUE})? If \code{FALSE} random counts or occurrences are generated using random draws from a suitably parameterised distribution, as stated in \code{countModel}.
##'
##' @return a matrix of simulated count or occurrence data, one row per gradient location, one column per species.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @keywords datagen
##'
##' @examples
##'
##' ## Poisson counts along a single gradient, Gaussian response
##'
##' x <- seq(from = 4, to = 6, length = 100)
##' opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
##' tol <- rep(0.25, 5)
##' h <- rep(20, 5)
##'
##' ## set.seed(1)
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson")
##' head(y)
##'
##' y <- coenocline(x, responseModel = "gaussian",
##'                 params = cbind(opt = opt, tol = tol, h = h),
##'                 countModel = "poisson",
##'                 expectation = TRUE)
##' matplot(x, y, type = "l", lty = "solid")
##'
##'
##' ## Poisson counts along two correlated gradients, Gaussian response
##' set.seed(1)
##' x1 <- runif(300, min = 4, max = 6)
##' opt1 <- seq(4, 6, length = 5)
##' tol1 <- rep(0.25, 5)
##' x2 <- runif(300, min = 2, max = 20)
##' opt2 <- seq(2, 20, length = 5)
##' tol2 <- rep(1, 5)
##' h <- rep(20, 5)
##'
##' set.seed(1)
##' params <- list(px = list(opt = opt1, tol = tol1, h = h),
##'                py = list(opt = opt2, tol = tol2))
##' y <- coenocline(cbind(x1, x2),
##'                 responseModel = "gaussian",
##'                 params = params,
##'                 extraParams = list(corr = 0.5),
##'                 countModel = "poisson")
##' head(y)
##'
##' ## Recreate Fig. 2 of Minchin (1987)
##' A0 <- c(5,4,7,5,9,8) * 10
##' m <- c(25,85,10,60,45,60)
##' r <- c(3,3,4,4,6,5) * 10
##' alpha <- c(0.1,1,2,4,1.5,1)
##' gamma <- c(0.1,1,2,4,0.5,4)
##' x <- 1:100
##' params <- list(m = m, A0 = A0, r = r, alpha = alpha, gamma = gamma)
##'
##' ## Expectations
##' set.seed(2)
##' y <- coenocline(x, responseModel = "beta",
##'                 params = params,
##'                 countModel = "poisson")
##' head(y)
##' matplot(y, type = "l", lty = "solid")
##'
##' y <- coenocline(x, responseModel = "beta",
##'                 params = params,
##'                 countModel = "poisson", expectation = TRUE)
##' matplot(y, type = "l", lty = "solid")
`coenocline` <- function(x,
                         responseModel = c("gaussian","beta"),
                         params,
                         extraParams = NULL,
                         countModel = c("poisson", "negbin", "bernoulli", "binary"),
                         countParams = NULL,
                         expectation = FALSE) {
    responseModel <- match.arg(responseModel)
    responseModel <- switch(responseModel,
                            gaussian = Gaussian,
                            beta     = Beta)
    countModel <- match.arg(countModel)
    countModel <- switch(countModel,
                         poisson = Poisson,
                         negbin  = NegBin,
                         binomial, binary = Bernoulli)

    ## x needs to be a vector, or for bivariate;
    ##   a list of 2 vectors, or a matrix of 2 columns
    ll <- is.list(x)
    mm <- is.matrix(x)

    if (ll | mm) { ## 2d
        if (ll) {
            stopifnot(length(x) == 2)
            X <- x[[1]]
            Y <- x[[2]]
        } else {
            stopifnot(NCOL(x) == 2)
            X <- x[, 1]
            Y <- x[, 2]
        }
        sim <- coenocline2d(X, Y,
                            responseModel = responseModel,
                            params = params, extraParams = extraParams,
                            countModel = countModel, countParams = countParams,
                            expectation = expectation)
    } else {
        sim <- coenocline1d(x,
                            responseModel = responseModel,
                            params = params, extraParams = extraParams,
                            countModel = countModel, countParams = countParams,
                            expectation = expectation)
    }
    sim
}

`coenocline1d` <- function(x, responseModel, params,
                           extraParams = NULL,
                           countModel, countParams = NULL,
                           expectation = FALSE) {
    n <- length(x)
    if (is.matrix(params)) {
        args <- list(x = x, params = params)
    } else if (is.list(params)) {
        if (is.list(params[[1]])) {
            params <- params[[1]]
        }
        args <- vector(mode = "list", length = length(params) + 1)
        args[[1]] <- x
        args[-1] <- params
        names(args) <- c("x", names(params))
    }
    ex <- do.call("expand", args)
    ## build arg list for responseModel
    if (is.null(extraParams)) {
        rargs <- list(x = ex[,1], px = as.list(data.frame(ex[, -1])))
    } else {
        le <- length(extraParams)
        rargs <- vector("list", length = 2 + le)
        names(rargs) <- c("x", "px", names(extraParams))
        rargs[["x"]] <- ex[,1]
        rargs[["px"]] <- as.list(data.frame(ex[, -1]))
        rargs[-(1:2)] <- extraParams
    }
    mu <- do.call(responseModel, rargs)
    if (expectation) {
        sim <- mu
    } else {
        if (is.null(countParams)) {
            cargs <- list(n = NROW(ex), mu = mu)
        } else {
            cargs <- vector("list", length = length(countParams) + 2)
            cargs[[1]] <- NROW(ex)
            cargs[[2]] <- mu
            cargs[-(1:2)] <- countParams
            names(cargs) <- c("n", "mu", names(countParams))
        }
        sim <- do.call(countModel, cargs)
    }
    sim <- matrix(sim, nrow = n)
    sim
}

`coenocline2d` <- function(x, y, responseModel, params,
                           extraParams = NULL,
                           countModel, countParams = NULL,
                           expectation = FALSE) {
    expandFun <- function(x, params) {
        args <- vector(mode = "list", length = length(params) + 1)
        args[[1]] <- x
        args[-1] <- params
        names(args) <- c("x", names(params))
        do.call("expand", args)
    }
    n1 <- length(x)
    n2 <- length(y)
    stopifnot(isTRUE(all.equal(n1, n2)))
    stopifnot(length(params) == 2L)
    exx <- expandFun(x, params[["px"]])
    exy <- expandFun(y, params[["py"]])
    if (is.null(extraParams)) {
        rargs <- list(x = exx[,1],
                      y = exy[,1],
                      px = as.list(data.frame(exx[, -1])),
                      py = as.list(data.frame(exy[, -1])))
    } else {
        le <- length(extraParams)
        rargs <- vector("list", length = 4 + le)
        names(rargs) <- c("x", "y", "px", "py", names(extraParams))
        rargs[["x"]] <- exy[,1]
        rargs[["y"]] <- exy[,1]
        rargs[["px"]] <- as.list(data.frame(exx[, -1]))
        rargs[["py"]] <- as.list(data.frame(exy[, -1]))
        rargs[-(1:4)] <- extraParams
    }
    mu <- do.call(responseModel, rargs)
    if (expectation) {
        sim <- mu
    } else {
        if (is.null(countParams)) {
            cargs <- list(n = NROW(exx), mu = mu)
        } else {
            cargs <- vector("list", length = length(countParams) + 2)
            cargs[[1]] <- NROW(exx)
            cargs[[2]] <- mu
            cargs[-(1:2)] <- countParams
            names(cargs) <- c("n", "mu", names(countParams))
        }
        sim <- do.call(countModel, cargs)
    }
    sim <- matrix(sim, nrow = n1)
    sim
}
