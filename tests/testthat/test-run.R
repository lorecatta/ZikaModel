test_that("run works", {
  set.seed(123)
  r1 <- run_deterministic_model()

  expect_type(r1$output, "double")
})

test_that("dt is used correctly to calculate time steps", {
  set.seed(123)
  r2 <- run_deterministic_model()

  mod <- r2$model
  results <- mod$transform_variables(r2$output)

  max_t <- r2$parameters$time_period / r2$parameters$DT
  tt <- seq(from = 1, to = max_t)
  expect_identical(tt * r2$parameters$DT, results$TIME)
})
