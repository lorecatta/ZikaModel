
devtools::load_all()

agec <- c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10)

death <- c(1e-10,
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

time <- 18200 # 50 years

set.seed(786808741)

model_run <- run_model(age,
                       death,
                       nn_links,
                       time = time)

save_plot(model_run$plot,
          "figures",
          "compartments_human",
          wdt = 17,
          hgt = 12)


# Scaling factor (between 0 and 1) for effect of seasonality on adult mosquitoes mortality. 1 = maximum effect of seasonality. 0 = no effect of seasonality. Default = 1.
# Kc_season Scaling factor (between 0 and 1) for effect of seasonality on mosquito larvae carrying capacity. Default = 1.
#  eip_season Scaling factor (between 0 and 1) for effect of seasonality on Extrinsic Incubation Period. Default = 1.
