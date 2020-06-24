test_that("format output works", {

  set.seed(123)
  r1 <- run_deterministic_model(time_period = 100)


  ### only human compartments


  ## summarise by compartment
  o1 <- format_output_H(r1)

  expect_s3_class(o1, "data.frame")
  expect_equal(sum(o1$compartment == "S"), 100)
  expect_equal(sum(o1$compartment == "I1"), 100)
  expect_equal(sum(o1$compartment == "R1"), 100)


  ## summarise by patch
  o2 <- format_output_H(r1, keep = "patch")
  expect_equal(dim(o2), c(6300, 4))

  ## summarise by vaccine status
  o3 <- format_output_H(r1, keep = "vaccine")
  expect_equal(dim(o3), c(600, 4))


  ### selected variables


  ## summarise by compartment
  o4 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"))
  expect_equal(dim(o4), c(400, 3))

  ## summarise by patch
  o5 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"), keep = "patch")
  expect_equal(dim(o5), c(8400, 4))

  ## summarise by vaccine status
  o6 <- format_output_H(r1, var_select = c("S", "R1", "inf_1", "MC_w"), keep = "vaccine")
  expect_equal(sum(o6$compartment == "S"), 200)
  expect_equal(sum(o6$compartment == "R1"), 200)
  expect_equal(sum(o6$compartment == "inf_1"), 200)
  expect_equal(sum(o6$compartment == "MC_w"), 200)


  #### mosquitoes


  ## summarise by compartment
  o7 <- format_output_M(r1)
  expect_equal(dim(o7), c(800, 3))

  ## summarise by patch
  o8 <- format_output_M(r1, keep = "patch")
  expect_equal(dim(o8), c(16800, 4))


  ### selected variables


  ## summarise by compartment
  o9 <- format_output_M(r1, var_select = c("Lwt", "Mwt_E2", "Mwb_S", "Delta"))
  expect_equal(sum(o9$compartment == "Lwt"), 100)
  expect_equal(sum(o9$compartment == "Mwt_E2"), 100)
  expect_equal(sum(o9$compartment == "Mwb_S"), 100)
  expect_equal(sum(o9$compartment == "Delta"), 100)


  ## summarise by patch
  o10 <- format_output_M(r1, var_select = c("Lwt", "Mwt_E2", "Mwb_S", "Delta"), keep = "patch")
  expect_equal(dim(o10), c(8400, 4))

})

test_that("format output mosquitoes proportions works", {

  r1 <- run_deterministic_model(time_period = 100)

  o1 <- format_output_Mprop(r1)
  o2 <- format_output_Mprop(r1, keep = "patch")

  expect_equal(dim(o1), c(100, 3))
  expect_equal(dim(o2), c(2100, 4))

})

test_that("summarising mosquito variables by vaccine status throws error", {

  r1 <- run_deterministic_model()

  expect_error(format_output_M(r1, keep = "vaccine"), "Can not summarise mosquito-related variables or compartments by vaccine status")

})
