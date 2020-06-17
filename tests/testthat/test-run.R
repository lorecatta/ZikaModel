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
  expect_identical(tt * r2$parameters$DT, results$time)
})

test_that("format output works", {

  set.seed(123)
  r1 <- run_deterministic_model()

  ## summarise by compartment
  o1 <- format_output_H(r1)

  n_comp <- length(levels(o1$compartment))
  expect_s3_class(o1, "data.frame")
  expect_equal(nrow(o1), n_comp * r1$parameters$time_period)

  ##summarise by patch
  o2 <- format_output_H(r1, keep = "patch")
  n_comp <- length(levels(o2$compartment))
  expect_equal(nrow(o2), n_comp * r1$parameters$time_period * r1$parameters$NP)

  ##summarise by vaccine status
  o3 <- format_output_H(r1, keep = "vaccine")
  n_comp <- length(levels(o3$compartment))
  expect_equal(nrow(o3), n_comp * r1$parameters$time_period * 2)

})
