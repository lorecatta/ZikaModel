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

test_that("format output works", {

  set.seed(123)
  r1 <- run_deterministic_model()

  ## summarise by compartment
  o1 <- format_output_H(r1)

  n_comp <- unique(o1$compartment)
  expect_type(o1, "data.frame")
  expect_equal(now(o1), n_comp * r2$parameters$time_period)

  ##summarise by patch
  browser()
  o2 <- format_output_H(r1, keep = "patch")

})
