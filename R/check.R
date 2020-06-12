
# -----------------------------------------------------------------------------

#' The function creates an equilibrium initialisation state to be used within
#' later model runs.
#'
#' @title Generate equilibrium state
#'
#' @param agec Vector of age group widths.
#'
#' @param death Numeric of mortality rates.
#'
#' @param model_parameter_list List of user-defined model parameters.
#'
#' @return Initial states of the model
#'
#' @export
equilibrium_init_create <- function(agec,
                                    death,
                                    model_parameter_list) {

  na <- as.integer(length(agec))  # number of age groups

  mpl <- model_parameter_list
  nn <- ZikaModel::nn_links # 8 nearest-neigbors to each patch.

  NP <- mpl$NP

  nn[nn > NP] <- NP

  surv <- ageb <- mean_surv <- mean_age <- lifeexpec <- corr_death <- c()

  surv[1] <- 1

  for (i in 2:(na+1)){

    #message("age group = ", i)
    surv[i] <- surv[i-1] * exp(-death[i-1] * agec[i-1])

  }

  ageb[1] <- 0

  for (i in 2:(na+1)){

    #message("age group = ", i)
    ageb[i] <- ageb[i-1] + agec[i-1]

  }

  for (i in 2:(na+1)){

    #message("age group = ", i)
    mean_surv[i-1] <- (surv[i-1] + surv[i]) * 0.5

  }

  for (i in 2:(na+1)){

    #message("age group = ", i)
    mean_age[i-1] <- (0.5 * ageb[i-1] + 0.5 * ageb[i])

  }

  for (i in 2:(na+1)){

    #message("age group = ", i)
    lifeexpec[i-1] <- (surv[i-1] - surv[i]) / death[i-1]

  }

  av_age_s <- agec * mean_surv

  lifespan <- sum(lifeexpec)

  corr_death[1] <- 1 / lifeexpec[1] - 1 / agec[1]

  for (i in 2:na){

    corr_death[i] <- (lifeexpec[i-1]/agec[i-1]) / lifeexpec[i] - 1 / agec[i]

  }

  YL <- mpl$YL
  DT <- mpl$DT

  deathrt <- corr_death / YL * DT

  N_human <- mpl$N_human

  N_eq <- c()

  N_eq[1:(NP-1)] <- N_human
  N_eq[NP] <- N_human / 10

  Nb <- N_eq / lifespan # number of births

  # initial compartment states

  other_prop_immune <- mpl$other_prop_immune
  other_foi <- mpl$other_foi

  aa <- bb <- c()

  for (i in 2:(na+1)){

    aa[i-1] <- (surv[i-1]-surv[i])/death[i-1]*(1-other_prop_immune*((1-exp(-other_foi*mean_age[i-1]))^1))

    bb[i-1] <- (surv[i-1]-surv[i])/death[i-1]*other_prop_immune*((1-exp(-other_foi*mean_age[i-1]))^1)

  }

  init_S <- array(0, c(na, 2, NP))

  init_S[1:na, 1, 1:NP] <- trunc(aa %o% Nb) # outer product

  init_S_sign <- signif(init_S, 6)

  init_I1 <- array(0, c(na, 2, NP))

  init_R1 <- array(0, c(na, 2, NP))

  init_R1[1:na, 1, 1:NP] <- trunc(bb %o% Nb)

  res <- list(nn = nn,
              na = na,
              surv = surv,
              mean_surv = mean_surv,
              mean_age = mean_age,
              av_age_s = av_age_s,
              deathrt = deathrt,
              lifeexpec = lifeexpec,
              lifespan = lifespan,
              N_eq = N_eq,
              Nb = Nb,
              init_S = init_S_sign,
              init_I1 = init_I1,
              init_R1 = init_R1)

  return(res)

}
