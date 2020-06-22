test_that("run works", {

  set.seed(123)
  r1 <- run_deterministic_model()

  expect_type(r1$output, "double")

})

test_that("dt is used correctly to calculate time steps", {

  set.seed(123)
  r1 <- run_deterministic_model()

  mod <- r1$model
  results <- mod$transform_variables(r1$output)

  max_t <- r1$parameters$time_period / r1$parameters$DT
  tt <- seq(from = 1, to = max_t)
  expect_identical(tt * r1$parameters$DT, results$time)

})

test_that("there are no infections when BG_FOI is 0", {

  r1 <- run_deterministic_model(BG_FOI = 0)

  o1 <- format_output_H(r1)
  expect_true(all(o1[o1$compartment == "I1", "y"] == 0))

  mod <- r1$model
  results <- mod$transform_variables(r1$output)
  t_max <- 365
  expect_equal(results$I1[t_max,,,], r1$parameters$init_I1)
  expect_equal(results$S[t_max,,,], r1$parameters$init_S)

})
