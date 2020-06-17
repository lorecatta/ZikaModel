


# -----------------------------------------------------------------------------

#' The function estimates the number of microcephaly cases due to ZIKV infection
#' during pregnancy.
#'
#' @title Estimates number of microcephaly cases caused by Zika infection
#'
#' @param x Zika_model_simulation object.
#'
#' @param var_select Vector of compartment names, e.g. \code{c("S", "R")}. In
#'   addition a number of additional variables can be requested. These include:
#' @export
calculate_microcases_ZIKV <- function(N_inf, pregnancy_risk_curve, birth_rates) {

  # pregnancy_risk_curve -


  # calculate the number of ZIKV infected foeti for each simulation time step

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

calculate_microcases_baseline <- function(N_tot, baseline_probM, birth_rates) {

  sweep(N_tot / 2, MARGIN = 2, birth_rates * baseline_probM, "*")

}

calculate_microcases <- function(N_inf, pregnancy_risk_curve, birth_rates, N_tot, baseline_probM) {

  a <- calculate_microcases_ZIKV(N_inf, pregnancy_risk_curve, birth_rates)

  b <- calculate_microcases_baseline(N_tot, baseline_probM, birth_rates)

  a + b

}

