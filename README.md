# zika-transmission-model
This repository contains the R package for running the Zika transmission model 
developed at Imperial College London using R odin.

The transmission model is a stochastic metapopulation model which includes the 
dynamics of the _Aedes_ mosquito vector and the age-stratified human host 
populations. The model simulates also the impact of seasonality and control 
strategies, such as the release of Wolbachia-infected mosquitoes and child 
vaccinationon, on virus dynamics .

For details of the original transmission model please see the 
[Ferguson et al. 2016 paper](https://science.sciencemag.org/content/353/6297/353) 
which is the article where the model is published.

## Installation
You need to first install the [odin](https://github.com/mrc-ide/odin) R package. 
The odin package allows you to write differential equations using a language 
similar to R, compile this into C code and find an approximate solution to the 
differential equations using a numerical method. 
Once odin is installed, you can install the `Zika_model` package by cloning this 
repository to your local machine, opening the `Zika_model.Rproj` file and load 
the `Zika_model` package using `devtools::load_all()`.

## Running the base model
To run the base version of the model, without seasonality and interventions, 
you can do the following:

```
# create a vector of human age groups 
age_init <- c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10)
  
# create a vector of human mortality rates 
deathrt <- c(1e-10, 
             1e-10, 
             1e-10, 
             0.00277068683332695, 
             0.0210680857689784,
             0.026724997685722,
             0.0525354529367476,
             0.0668013582441452,
             0.119271483740379,
             0.279105747097929,
             0.390197266957464)
             
# provide the length of time (in days) that you want to run the model for
time_frame <- 364 * 50
  
# run the model
model_run <- run_model(age = age_init,
                       death = deathrt,
                       nn_links,
                       time = time_frame)
```
  
You can use `save_plot` to save a png figure of the plot of human compartments 
against time from the model run.

```
save_plot(plot_obj = model_run$plot, out_pth = "figures", out_fl_nm = "compartments_human", wdt = 17, hgt = 12)
```

## Seasonality 
The model allows to account for the effect of seasonal variations in climatic 
variables (e.g. temperature and precipitation) on Zika transmission dynamics. 
In the model, seasonality affects: 

1. adult mosquitoes mortality;
2. mosquito larvae carrying capacity;
3. Extrinsic Incubation Period.

At the moment the effect of seasonality is implemented as "all or nothing" 
(e.g. it can either be switched completely on or off). 
Intermediate effects are possible but they require a (simple) change to the code.

A model which includes the effect of seasonality can be implemnted as:

```
seasonal_model_run <- run_model(age = age_init,
                                death = deathrt,
                                nn_links,
                                time = time_period,
                                season = TRUE)
```
