
#' @noRd
calculate_mc_probability_given_ZIKV_infection_pregnancy <- function() {

  ## Script for Lorenzo to look at microcephaly risk curves

  ## devtools::install_github("jameshay218/zikaInfer")

  ## Reads in all files ending "multivariate_chain" from specified working dir
  chains <- zikaInfer::load_mcmc_chains(location = ".", asList = FALSE, unfixed = FALSE)

  ## Gamma function used
  dput(microceph_v1)

  chain_means <- colMeans(chains)
  gamma_mean <- as.numeric(chain_means["mean"])
  gamma_var <- as.numeric(chain_means["var"])
  gamma_c <- as.numeric(chain_means["c"])

  pars <- c(mean = gamma_mean,
            var = gamma_var,
            c = gamma_c)

  # probability that a foetus developed microcephaly
  # given that the mother was infected in a particular week during pregnancy
  m_risk_probs <- microceph_v1(pars = pars)

  data.frame(day = 0:279, prob = m_risk_probs)

}


#' @noRd
calculate_brazil_birth_rates <- function(br_brazil_all_years) {

  br_brazil <- br_brazil_all_years[, "X2010...2015"] / 1000

  br_brazil_2 <- c(0, br_brazil)

  br_brazil_av <- c()

  my_seq <- seq(1, length(br_brazil_2), by = 2)

  for (i in seq_along(my_seq)) {

    index <- my_seq[i]

    out <- (br_brazil_2[index] + br_brazil_2[index + 1]) / 2

    br_brazil_av[i] <- out

  }

  br_brazil_age <- rep(0, na)

  br_brazil_age[1] <- 0
  br_brazil_age[2] <- 0
  br_brazil_age[3] <- br_brazil_av[1]
  br_brazil_age[4] <- br_brazil_av[2]
  br_brazil_age[5] <- br_brazil_av[3]
  br_brazil_age[6] <- br_brazil_av[4]
  br_brazil_age[7] <- 0
  br_brazil_age[8] <- 0
  br_brazil_age[9] <- 0
  br_brazil_age[10] <- 0
  br_brazil_age[11] <- 0

  # br_brazil_age_2 <- br_brazil_age / YL * DT

  br_brazil_age

}

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

