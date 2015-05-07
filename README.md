coenocliner
===========

#### Released version
[![CRAN version](http://www.r-pkg.org/badges/version/coenocliner)](http://cran.rstudio.com/web/packages/coenocliner/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/coenocliner)](http://cran.rstudio.com/web/packages/coenocliner/index.html)

#### Build status
[![Build Status](https://travis-ci.org/gavinsimpson/coenocliner.svg?branch=master)](https://travis-ci.org/gavinsimpson/coenocliner)  [![Build status](https://ci.appveyor.com/api/projects/status/hc8dbxrim2nj3c1i/branch/master)](https://ci.appveyor.com/project/gavinsimpson/coenocliner/branch/master)

## What is coenocliner?

An R package to simulate species abundances (counts) along gradients

One of the key ways quantitative ecologists attempt to understand 
the properties and behaviour of the methods they use or dream up is 
through the use of simulated data. There are a number of computer 
programmes for simulating ecological data along gradients, such as 
Peter Minchin's COMPAS, but none (that I am aware of) that are 
available for R on CRAN. Dave Robert's 
[coenoflex](http://cran.r-project.org/package=coenoflex) package for 
R would be a useful alternative but currently is archived on CRAN 
because of some problems in the Fortran code underlying the package.

Rather than have to reinvent the wheel each time I wanted to simulate 
some new data for a paper or to work on a new approach, I decided to 
start my own R package to contain a range of simulators encapsulating 
different response models, numbers of gradients, etc.

At the moment, coenocliner is limited in what it can do practically. 
There is a single response model, the Gaussian response, which is a 
symmetric model of the parameters; the optimum, tolerance and height 
of the response curve. Count data can be generated from this model 
from either a Poisson or negative binomial distribution, using the 
parameterised Gaussian response as the expectation or mean of the 
distribution.

Additional response models include:

 1. The generalised beta response function

A further feature of **coenocliner** that I hope to develop is to 
include simulation  wrapper functions that replicate the simulation 
methods used in research papers. A working example is `simJamil`, 
which produces simlations from a Gaussian logit response following 
the scheme described in Jamil & ter Braak (2013).

## Development

I would like to see coenocliner be as inclusive as possible; if you 
have code to simulate ecological species or community data that is 
just sitting around, consider adding it to coenocliner. In the 
meantime, I'm happy just having something tangible for my own use 
without having to remember the expressions for some of the response 
models.

Currently coenocliner is licensed under the GPL v2, but I'm happy to 
reconsider this if you want to contribute code under a more permissive 
licence.

## Installation

No binary packages are currently available for coenocliner. If you 
have the correct development tools you can compile the package 
yourself after downloading the source code from github. Once I work 
out how to link git with svn I'll start a project on 
[R-forge](http://r-forge.r-project.org) which will host binary 
packages of coenocliner.

If you use Hadley Wickham's **devtools** package then you 
can install coenocliner directly from github using functions that 
devtools provides. To do this, install **devtools** from CRAN via

    install.packages("devtools")

then run

    devtools::install_github("gavinsimpson/coenocliner")

### References

Jamil and ter Braak (2013) Generalized linear mixed models can 
detect unimodal species-environment relationships. *PeerJ* **1:e95**;
 [DOI 10.7717/peerj.95](http://doi.org/10.7717/peerj.95)
