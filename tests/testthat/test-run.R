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


  ### only human compartments


  ## summarise by compartment
  o1 <- format_output_H(r1)

  n_comp <- length(levels(o1$compartment))
  expect_s3_class(o1, "data.frame")
  expect_true(nrow(o1) == n_comp * r1$parameters$time_period)

  ## summarise by patch
  o2 <- format_output_H(r1, keep = "patch")
  n_comp <- length(levels(o2$compartment))
  expect_equal(nrow(o2), n_comp * r1$parameters$time_period * r1$parameters$NP)

  ## summarise by vaccine status
  o3 <- format_output_H(r1, keep = "vaccine")
  n_comp <- length(levels(o3$compartment))
  expect_equal(nrow(o3), n_comp * r1$parameters$time_period * 2)


  ### selected variables


  ## summarise by compartment
  o4 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"))
  n_comp <- length(levels(o4$compartment))
  expect_equal(nrow(o4), n_comp * r1$parameters$time_period)

  ## summarise by patch
  o5 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"), keep = "patch")
  n_comp <- length(levels(o5$compartment))
  expect_equal(nrow(o5), n_comp * r1$parameters$time_period * r1$parameters$NP)

  ## summarise by vaccine status
  o6 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"), keep = "vaccine")
  n_comp <- length(levels(o6$compartment))
  expect_equal(nrow(o6), n_comp * r1$parameters$time_period * 2)


  #### mosquitoes


  ## summarise by compartment
  o7 <- format_output_M(r1)
  n_comp <- length(levels(o7$compartment))
  expect_equal(nrow(o7), n_comp * r1$parameters$time_period)

  ## summarise by patch
  o8 <- format_output_M(r1, keep = "patch")
  n_comp <- length(levels(o8$compartment))
  expect_equal(nrow(o8), n_comp * r1$parameters$time_period * r1$parameters$NP)


  ### selected variables


  ## summarise by compartment
  o9 <- format_output_M(r1, var_select = c("Lwt", "Mwt_E2", "Mwb_S", "Delta"))
  n_comp <- length(levels(o9$compartment))
  expect_equal(nrow(o9), n_comp * r1$parameters$time_period)

  ## summarise by patch
  o10 <- format_output_M(r1, var_select = c("Lwt", "Mwt_E2", "Mwb_S", "Delta"), keep = "patch")
  n_comp <- length(levels(o10$compartment))
  expect_equal(nrow(o10), n_comp * r1$parameters$time_period * r1$parameters$NP)

})
