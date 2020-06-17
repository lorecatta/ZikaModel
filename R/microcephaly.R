
# -----------------------------------------------------------------------------

#' The function estimates the number of microcephaly cases due to ZIKV infection
#' during pregnancy.
#'
#' @title Estimates number of microcephaly cases caused by Zika infection
#'
#' @param model Zika_model_simulation object.
#'
#' @param birth_rates age specific birth rates for Brazil.
#'
#' @export
calculate_microcases_ZIKV <- function(model, birth_rates) {

  N_inf <- unpack_odin(model, "inf_1")

  pregnancy_risk_curve <- ZikaModel::mc_prob_ZIKV_pregn

  array_dim <- dim(N_inf)

  probM <- pregnancy_risk_curve$prob

  n_time_steps <- array_dim[1]

  probM_size <- length(probM)

  all_probs <- array(0, array_dim)

  for (i in seq(1, n_time_steps, 1)) {

    tmp <- 0

    testDay <- i

    minDay <- testDay - probM_size + 1

    if(minDay < 1) minDay <- 1

    #message("testDay = ", testDay)
    #message("minDay = ", minDay)

    for(j in seq(minDay, testDay, by = 1)) {

      #message("j = ", j)

      index <- j - (testDay - probM_size)
      #message("index = ", index)

      if(index == 0) stop("j = ", j)

      # tmp <- tmp + (N_inf[j,,,]/2) * (probM[index] + bp - probM[index]*bp)
      tmp <- tmp + (N_inf[j,,,] / 2) * probM[index]
    }

    #message("tmp = ", tmp)

    all_probs[i,,,] <- tmp #1 - (1 - tmp) * (1 - bp)

  }

  sweep(all_probs, MARGIN = 2, birth_rates, "*")

}


# -----------------------------------------------------------------------------

#' The function estimates the baseline number of microcephaly cases
#' (not due to ZIKV infection during pregnancy).
#'
#' @title Estimates baseline number of microcephaly cases
#'
#' @param model Zika_model_simulation object.
#'
#' @param birth_rates age specific birth rates for Brazil.
#'
#' @export
calculate_microcases_baseline <- function(model, birth_rates) {

  N_tot <- unpack_odin(model, "Ntotal")

  baseline_probM <- model$parameters$mc_baseline

  sweep(N_tot / 2, MARGIN = 2, birth_rates * baseline_probM, "*")

}


# -----------------------------------------------------------------------------

#' The function estimates the total number of microcephaly cases
#' (baseline and due to ZIKV infection during pregnancy).
#'
#' @title Estimates total number of microcephaly cases
#'
#' @param model Zika_model_simulation object.
#'
#' @export
calculate_microcases <- function(model) {

  birth_rates <- ZikaModel::br_brazil_age

  birth_rates_2 <- birth_rates / model$parameters$YL * model$parameters$DT

  a <- calculate_microcases_ZIKV(model, birth_rates_2)

  b <- calculate_microcases_baseline(model, birth_rates_2)

  a + b

}
