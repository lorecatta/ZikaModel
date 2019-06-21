
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

1.  First install `devtools`, if you don't already have it.

``` r
install.packages("devtools")
library(devtools)
```

1.  Then, in a fresh R session, install the `ZikaModel` package.

``` r
devtools::install_github("mrc-ide/zika-transmission-model")
#> Skipping install of 'ZikaModel' from a github remote, the SHA1 (29716090) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

1.  Load and attach it

``` r
library(ZikaModel)
```

Running the base model
----------------------

To run the base version of the model, without seasonality and interventions, you can do the following:

``` r
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
  
set.seed(1)

# run the model
model_run <- ZikaModel::run_model(agec = age_init,
                                  death = deathrt,
                                  nn_links,
                                  time = time_frame)
#> Equations use index variables i on the rhs outside of an index.
#> The behaviour of this has changed since odin 0.1.3 - see
#> https://github.com/mrc-ide/odin/issues/136 for details.
#> To silence this note, set option `odin.no_check_naked_index` to TRUE
#> This note will disappear in a version after odin 1.0.0
#>  vacc_noncov[1:vnc_row, 1] <- (if ((step >= YL*vacc_child_starttime) && # (line 125)
#>                                    (step < YL*vacc_child_stoptime) && # (line 126)
#>                                    ((i == vacc_child_age + 1))) # (line 127)
#>    (1 - vacc_child_coverage) else 1) * # (line 128)
#>    (if ((step == vacc_cu_rndtime) && # (line 129)
#>         ((i >= vacc_cu_minage + 1)) && # (line 130)
#>         ((i <= vacc_cu_maxage + 1))) (1 - vacc_cu_coverage) else 1) # (line 131)
#>  vacc_noncov[1:vnc_row, 2] <- 1 # (line 133)
#> Unused equations: age_inf1, AGE_REC, dim_mean_age, dis_scale, dis_scaleM, dis_scaleW, disease_patch, disease_patch_cum, mean_age, mean_age_inf1_nv, mean_age_inf1_v, NTnv, Ntotal_np, Ntotal_nv, Ntotal_v, NTv, overall_mean_age_inf1, p_age_dist1, p_age_inf1, p_sum_inf1, prop_wb, PropDiseaseReported, propTransGlobal_bigpatch, R1nv, Snv, sum_age_inf1, sum_inf1
#>  age_inf1[,,] <- sinf1[i,j,k] * mean_age[i] # (line 594)
#>  AGE_REC <- user() # (line 437)
#>  dim(mean_age) <- na # (line 831)
#>  dis_scale <- PropDiseaseReported * 10000 # (line 453)
#>  dis_scaleM <- dis_scale * YL / 30          # monthly rates # (line 455)
#>  dis_scaleW <- dis_scale                    # weekly rates # (line 454)
#>  disease_patch[] <- disease_patch_cum[na,i] # (line 481)
#>  disease_patch_cum[1,1:NP] <- disease_age_inc[i,1,j] + disease_age_inc[i,2,j] # (line 476)
#>  disease_patch_cum[2:na,1:NP] <- disease_patch_cum[i-1,j] + # (line 477)
#>                                  disease_age_inc[i,1,j] + # (line 478)
#>                                  disease_age_inc[i,2,j] # (line 479)
#>  mean_age[] <- user() # (line 593)
#>  mean_age_inf1_nv[] <- sum_age_inf1[na,1,i] / (1e-20 + sum_inf1[na,1,i]) # (line 602)
#>  mean_age_inf1_v[] <- sum_age_inf1[na,2,i] / (1e-20 + sum_inf1[na,2,i]) # (line 604)
#>  NTnv <- sum(Ntotal_nv[]) + 1e-20 # (line 334)
#>  Ntotal_np[1:na,1:2,1] <- Ntotal[i,j,k] # (line 325)
#>  Ntotal_np[1:na,1:2,2:NP] <- Ntotal_np[i,j,k-1] + Ntotal[i,j,k] # (line 326)
#>  Ntotal_nv[] <- Ntotal_np[i,1,NP] # (line 332)
#>  Ntotal_v[] <- Ntotal_np[i,2,NP] # (line 333)
#>  NTv <- sum(Ntotal_v[]) + 1e-20 # (line 335)
#>  overall_mean_age_inf1 <- sum(age_inf1[,,]) / (1e-20 + sum(sinf1[,,])) # (line 606)
#>  p_age_dist1[] <- p_age_inf1[i] / (1e-20 + sum(p_age_inf1[])) # (line 612)
#>  p_age_inf1[] <- p_sum_inf1[i,NP-1] # (line 611)
#>  p_sum_inf1[1:na,1] <- sinf1[i,1,j] + sinf1[i,2,j] # (line 608)
#>  p_sum_inf1[1:na,2:(NP-1)] <- p_sum_inf1[i,j-1] + sinf1[i,1,j] + sinf1[i,2,j] # (line 609)
#>  prop_wb[] <- Mwb_tot[i] / (M_tot[i] + 1e-10) # (line 232)
#>  PropDiseaseReported <- user() # (line 452)
#>  propTransGlobal_bigpatch <- propTransGlobal / 10 # (line 411)
#>  R1nv[,] <- R1[i,1,j] / Ntotal[i,1,j] # (line 338)
#>  Snv[,] <- S[i,1,j] / Ntotal[i,1,j] # (line 337)
#>  sum_age_inf1[1,1:2,1:NP] <- age_inf1[i,j,k] # (line 599)
#>  sum_age_inf1[2:na,1:2,1:NP] <- age_inf1[i,j,k] + sum_age_inf1[i-1,j,k] # (line 600)
#>  sum_inf1[1,1:2,1:NP] <- sinf1[i,j,k] # (line 596)
#>  sum_inf1[2:na,1:2,1:NP] <- sinf1[i,j,k] + sum_inf1[i-1,j,k] # (line 597)
```

You can use `save_plot` to save a png figure of the plot of human compartments against time from the model run.

``` r
ZikaModel::save_plot(plot_obj = model_run$plot, 
                     out_pth = "figures", 
                     out_fl_nm = "compartments_human", 
                     wdt = 17, 
                     hgt = 12)
```

At this point it might be useful to inspect some diagnostics to check that the model is actually doing what we want. The function `post_processing` reshapes the model outputs in a more convenient format and saves plots of selected diagnostics.

``` r
ZikaModel::post_processing(model_run$dat, "figures")
```

Seasonality
-----------

The model allows to account for the effect of seasonal variations in climatic variables (e.g. temperature and precipitation) on Zika transmission dynamics. In the model, seasonality affects:

-   adult mosquitoes mortality;
-   mosquito larvae carrying capacity;
-   Extrinsic Incubation Period.

At the moment the effect of seasonality is implemented as "all or nothing" (e.g. it can either be switched completely on or off). Intermediate effects are possible but they require a (simple) change to the code.

A model which includes the effect of seasonality can be implemented as:

``` r
set.seed(10)

seasonal_model_run <- ZikaModel::run_model(agec = age_init,
                                           death = deathrt,
                                           nn_links,
                                           time = time_frame,
                                           season = TRUE)
#> Equations use index variables i on the rhs outside of an index.
#> The behaviour of this has changed since odin 0.1.3 - see
#> https://github.com/mrc-ide/odin/issues/136 for details.
#> To silence this note, set option `odin.no_check_naked_index` to TRUE
#> This note will disappear in a version after odin 1.0.0
#>  vacc_noncov[1:vnc_row, 1] <- (if ((step >= YL*vacc_child_starttime) && # (line 125)
#>                                    (step < YL*vacc_child_stoptime) && # (line 126)
#>                                    ((i == vacc_child_age + 1))) # (line 127)
#>    (1 - vacc_child_coverage) else 1) * # (line 128)
#>    (if ((step == vacc_cu_rndtime) && # (line 129)
#>         ((i >= vacc_cu_minage + 1)) && # (line 130)
#>         ((i <= vacc_cu_maxage + 1))) (1 - vacc_cu_coverage) else 1) # (line 131)
#>  vacc_noncov[1:vnc_row, 2] <- 1 # (line 133)
#> Unused equations: age_inf1, AGE_REC, dim_mean_age, dis_scale, dis_scaleM, dis_scaleW, disease_patch, disease_patch_cum, mean_age, mean_age_inf1_nv, mean_age_inf1_v, NTnv, Ntotal_np, Ntotal_nv, Ntotal_v, NTv, overall_mean_age_inf1, p_age_dist1, p_age_inf1, p_sum_inf1, prop_wb, PropDiseaseReported, propTransGlobal_bigpatch, R1nv, Snv, sum_age_inf1, sum_inf1
#>  age_inf1[,,] <- sinf1[i,j,k] * mean_age[i] # (line 594)
#>  AGE_REC <- user() # (line 437)
#>  dim(mean_age) <- na # (line 831)
#>  dis_scale <- PropDiseaseReported * 10000 # (line 453)
#>  dis_scaleM <- dis_scale * YL / 30          # monthly rates # (line 455)
#>  dis_scaleW <- dis_scale                    # weekly rates # (line 454)
#>  disease_patch[] <- disease_patch_cum[na,i] # (line 481)
#>  disease_patch_cum[1,1:NP] <- disease_age_inc[i,1,j] + disease_age_inc[i,2,j] # (line 476)
#>  disease_patch_cum[2:na,1:NP] <- disease_patch_cum[i-1,j] + # (line 477)
#>                                  disease_age_inc[i,1,j] + # (line 478)
#>                                  disease_age_inc[i,2,j] # (line 479)
#>  mean_age[] <- user() # (line 593)
#>  mean_age_inf1_nv[] <- sum_age_inf1[na,1,i] / (1e-20 + sum_inf1[na,1,i]) # (line 602)
#>  mean_age_inf1_v[] <- sum_age_inf1[na,2,i] / (1e-20 + sum_inf1[na,2,i]) # (line 604)
#>  NTnv <- sum(Ntotal_nv[]) + 1e-20 # (line 334)
#>  Ntotal_np[1:na,1:2,1] <- Ntotal[i,j,k] # (line 325)
#>  Ntotal_np[1:na,1:2,2:NP] <- Ntotal_np[i,j,k-1] + Ntotal[i,j,k] # (line 326)
#>  Ntotal_nv[] <- Ntotal_np[i,1,NP] # (line 332)
#>  Ntotal_v[] <- Ntotal_np[i,2,NP] # (line 333)
#>  NTv <- sum(Ntotal_v[]) + 1e-20 # (line 335)
#>  overall_mean_age_inf1 <- sum(age_inf1[,,]) / (1e-20 + sum(sinf1[,,])) # (line 606)
#>  p_age_dist1[] <- p_age_inf1[i] / (1e-20 + sum(p_age_inf1[])) # (line 612)
#>  p_age_inf1[] <- p_sum_inf1[i,NP-1] # (line 611)
#>  p_sum_inf1[1:na,1] <- sinf1[i,1,j] + sinf1[i,2,j] # (line 608)
#>  p_sum_inf1[1:na,2:(NP-1)] <- p_sum_inf1[i,j-1] + sinf1[i,1,j] + sinf1[i,2,j] # (line 609)
#>  prop_wb[] <- Mwb_tot[i] / (M_tot[i] + 1e-10) # (line 232)
#>  PropDiseaseReported <- user() # (line 452)
#>  propTransGlobal_bigpatch <- propTransGlobal / 10 # (line 411)
#>  R1nv[,] <- R1[i,1,j] / Ntotal[i,1,j] # (line 338)
#>  Snv[,] <- S[i,1,j] / Ntotal[i,1,j] # (line 337)
#>  sum_age_inf1[1,1:2,1:NP] <- age_inf1[i,j,k] # (line 599)
#>  sum_age_inf1[2:na,1:2,1:NP] <- age_inf1[i,j,k] + sum_age_inf1[i-1,j,k] # (line 600)
#>  sum_inf1[1,1:2,1:NP] <- sinf1[i,j,k] # (line 596)
#>  sum_inf1[2:na,1:2,1:NP] <- sinf1[i,j,k] + sum_inf1[i-1,j,k] # (line 597)

ZikaModel::save_plot(plot_obj = seasonal_model_run$plot, 
                     out_pth = "figures/seasonality", 
                     out_fl_nm = "compartments_human", 
                     wdt = 17, 
                     hgt = 12)

ZikaModel::post_processing(seasonal_model_run$dat, "figures/seasonality")
```
