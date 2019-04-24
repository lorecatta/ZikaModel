# stochastic-zika-model
This repository contains the R package for running the Zika transmission model developed at Imperial College London using R odin.

The transmission model is a stochastic metapopulation model which includes the dynamics of the _Aedes_ mosquito vector and the human host populations. The model simulates also the impact of control strategies such as the release of Wolbachia-infected mosquitoes and child vaccination.

For details of the original transmission model please see the [Ferguson et al. 2016 paper](https://science.sciencemag.org/content/353/6297/353) 
which is the article where the model is published.

## Running the model
To run the model you can do the following:

```
# load the package (or have it built and reloaded as described above)
library(stochZika)

# create a vector of age groups 
age_init <- c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10)
  
# create a vector of mortality rates 
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
  
You can save a plot of the human compartments

```
save_plot(model_run$plot, "figures", "compartments_human", wdt = 17, hgt = 12)
```

The model allows to account for the effect of seasonal variations in climatic 
variables (e.g. temperature and precipitation) on Zika transmission dynamics. 
In the model, seasonality affects: 

1. adult mosquitoes mortality;
2. mosquito larvae carrying capacity;
3. Extrinsic Incubation Period.

At the moment the effect of seasonality is implemnted as "all or nothing" 
(e.g. it can either be switched completely on or off). 
Intermediate effects are possible but they require a (simple) change in the code.

A model which includes the effect of seasonality can be implemnted as:

```
seasonal_model_run <- run_model(age = age_init,
                                death = deathrt,
                                nn_links,
                                time = time_period,
                                season = TRUE)
```
