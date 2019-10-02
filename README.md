
<!-- README.md is generated from README.Rmd. Please edit that file -->
ZikaModel
=========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/mrc-ide/zika-transmission-model.svg?branch=master)](https://travis-ci.org/mrc-ide/zika-transmission-model) <!-- badges: end -->

`ZikaModel` is an R package for running the Zika transmission model developed at Imperial College London using R odin.

The transmission model is a metapopulation model which includes the dynamics of the *Aedes* mosquito vector and the age-stratified human host populations. The model which has a stochastic and a deterministic version, simulates also the impact of seasonality and control strategies, such as the release of Wolbachia-infected mosquitoes and child vaccinationon, on virus dynamics .

For details of the original transmission model please see the [Ferguson et al. 2016 paper](https://science.sciencemag.org/content/353/6297/353) which is the article where the model is published.

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

Check this vignette on how to run the [deterministic]() base version of the model.
