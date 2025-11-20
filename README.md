# coenocliner

<!-- badges: start -->
[![CRAN version](http://www.r-pkg.org/badges/version/coenocliner)](https://cran.r-project.org/package=coenocliner) [![](https://cranlogs.r-pkg.org/badges/grand-total/coenocliner)](https://cran.r-project.org/package=coenocliner)

[![R-CMD-check](https://github.com/gavinsimpson/coenocliner/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gavinsimpson/coenocliner/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/gavinsimpson/coenocliner/graph/badge.svg)](https://app.codecov.io/gh/gavinsimpson/coenocliner)
<!-- badges: end -->

## What is coenocliner?

An R package to simulate species abundances (counts) along gradients

One of the key ways quantitative ecologists attempt to understand 
the properties and behaviour of the methods they use or dream up is 
through the use of simulated data. There are a number of computer 
programmes for simulating ecological data along gradients, such as 
Peter Minchin's COMPAS, but none (that I am aware of) that are 
available for R on CRAN. Dave Robert's 
[coenoflex](https://cran.r-project.org/package=coenoflex) package for 
R is a useful alternative but at the time I started writing this, coenoflex was 
archived on CRAN (it is now back on CRAN).

Rather than have to reinvent the wheel each time I wanted to simulate 
some new data for a paper or to work on a new approach, I decided to 
start my own R package to contain a range of simulators encapsulating 
different response models, numbers of gradients, etc.

At the moment, **coenocliner** is limited in what it can do practically. 
There are two response models, the Gaussian response, which is a 
symmetric model of the parameters (the optimum, tolerance and height 
of the response curve) and the generalized beta response model.
Count data can be generated from this model from several distributions, including the Poisson or negative binomial distributions, using the 
parameterised response model as the expectation or mean of the 
distribution. The complete list of distributions supported are:

* Poisson,
* Negative binomial,
* Bernoulli,
* Binomial,
* Beta binomial,
* Zero-inflated Poisson,
* Zero-inflated negative binomial,
* Zero-inflated binomial,
* Zero-inflated beta binomial.

A further feature of **coenocliner** that I hope to develop is to 
include simulation wrapper functions that replicate the simulation 
methods used in research papers. A working example is `simJamil`, 
which produces simulations from a Gaussian logit response following 
the scheme described in Jamil & ter Braak (2013).

## Development

I would like to see **coenocliner** be as inclusive as possible; if you 
have code to simulate ecological species or community data that is 
just sitting around, consider adding it to **coenocliner**. In the 
meantime, I'm happy just having something tangible for my own use 
without having to remember the expressions for some of the response 
models.

Currently **coenocliner** is licensed under the GPL v2, but I'm happy to 
reconsider this if you want to contribute code under a more permissive 
licence.

## Installation

You can install the released version directly from CRAN using

```r
install.packages("coenocliner")
```

If you use the **remotes** package then you  can install **coenocliner**
directly from GitHub using functions that **remotes** provides. To do this, 
install **remotes** from CRAN via

    install.packages("remotes")

then run

    remotes::install_github("gavinsimpson/coenocliner")

### References

Jamil and ter Braak (2013) Generalized linear mixed models can 
detect unimodal species-environment relationships. *PeerJ* **1:e95**;
 [DOI 10.7717/peerj.95](https://doi.org/10.7717/peerj.95)
