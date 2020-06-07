test_that("run works", {

  set.seed(123)

  r1 <- run_deterministic_model()

})

test_that("dt is used correctly to calculate time steps", {

  age_init <- c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10)

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

  time_frame <- 364 * 10

  tt <- seq(0, time_frame, 1)

  mpl <- create_model_param_list(season = FALSE)

  dt <- mpl$DT

  model_outputs <- run_model(agec = age_init,
                             death = deathrt,
                             nn_links,
                             amplitudes_phases,
                             time = time_frame)

  expect_identical(tt * dt, model_outputs$TIME)

})
