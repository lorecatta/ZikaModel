
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

set.seed(10)

seasonal_model_run <- run_model(agec,
                                death,
                                nn_links,
                                time = time,
                                season = TRUE)

save_plot(seasonal_model_run$plot,
          "figures/seasonality",
          "compartments_human",
          wdt = 17,
          hgt = 12)

post_processing(seasonal_model_run$dat, "figures/seasonality")
