
<!-- README.md is generated from README.Rmd. Please edit that file -->
ZikaModel
=========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/mrc-ide/ZikaModel.svg?branch=master)](https://travis-ci.org/mrc-ide/ZikaModel) [![Codecov test coverage](https://codecov.io/gh/mrc-ide/ZikaModel/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/ZikaModel?branch=master) [![R build status](https://github.com/mrc-ide/ZikaModel/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/ZikaModel/actions) <!-- badges: end -->

`ZikaModel` is an R package for running the Zika transmission model developed at Imperial College London using R odin.

The transmission model is a metapopulation model which includes the dynamics of the *Aedes aegypti* mosquito vector and the age-stratified human host populations. The model which has a stochastic and a deterministic version, simulates also the effect of seasonality and the impact of control strategies, such as the release of Wolbachia-infected mosquitoes and child vaccination.

For details of the original transmission model please see the original [article](https://science.sciencemag.org/content/353/6297/353) where the model is published.

Installation
------------

You need to first install the [odin](https://github.com/mrc-ide/odin) R package.

Once odin is installed, you can install the `ZikaModel` package with the following steps:

-   First install `devtools`, if you don't already have it

``` r
install.packages("devtools")
library(devtools)
```

-   Then, in a fresh R session, install the `ZikaModel` package

``` r
devtools::install_github("mrc-ide/ZikaModel")
```

-   Load and attach it

``` r
library(ZikaModel)
```

Running the base model
----------------------

Check out this vignette on how to run the [deterministic](https://mrc-ide.github.io/ZikaModel/articles/deterministic_base.html) base version of the model.
