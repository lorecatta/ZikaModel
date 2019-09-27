
#------------------------------------------------------------------------------

# equilibrium_init_create

#' \code{equilibrium_init_create} creates an equilibrium initialisation
#' state to be used within later model runs.
#'
#' @param agec Vector of age group widths.
#'
#' @param death Numeric of mortality rates.
#'
#' @param nn_links 8 nearest-neigbors to each patch.
#'
#' @param model_parameter_list list of user-defined model parameters.

#' @importFrom stats rnorm
#'
#' @export


equilibrium_init_create <- function(agec, death, nn_links, model_parameter_list){

  mpl <- model_parameter_list
  nn <- nn_links

  NP <- mpl$NP

  nn[nn > NP] <- NP

  na <- as.integer(length(agec))  # number of age groups

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

  Mwt_cv <- mpl$Mwt_cv
  Mwt_mean <- mpl$Mwt_mean

  ln_sd <- sqrt(log(1 + Mwt_cv * Mwt_cv))
  ln_mean <- log(Mwt_mean) - (ln_sd * ln_sd) / 2

  init_Mwt_base <- exp(rnorm(NP, ln_mean, ln_sd))
  # init Mwt_base[1..NP]=exp(normal(ln_mean,ln_sd,seed+i))



  # -----------------------------------------------------------------------------
  #
  # Climate seasonality
  #
  # -----------------------------------------------------------------------------



  # Effect of spatial variation in seasonality on dengue transmission.

  # Patches have different seasonality of dengue transmission, depending on
  # where they are located with respect to the equator.

  # In the model seasonal climate variation affects:
  # 1 - Larval carrying capacity
  # 2 - Vector mortality
  # 3 - Extrinsic Incubation Period

  # Delta - vector mortality
  # Kc - carrying capacity
  # eip - extrinsic incubation period

  season_amp <- season_phase <- Wb_introtime <- c()

  season_phase[1:8] <- 0.5 * YL
  season_phase[9:12] <- 0.25 * YL
  season_phase[13:20] <- 0
  season_phase[21] <- 0.25 * YL

  season_amp[1:8] <- 1
  season_amp[13:20] <- 1
  season_amp[9:12] <- 0.33
  season_amp[21] <- 0.33

  Wb_starttime <- mpl$Wb_starttime

  for (i in seq_len(NP - 1)){

    Wb_introtime[i] <- Wb_starttime + trunc(6 * (i - 1) / NP) * 0.5

  }

  Wb_introtime[NP] <- 1000

  vacc_cu_time <- mpl$vacc_cu_time
  age_per <- mpl$age_per

  vacc_cu_rndtime <- trunc(vacc_cu_time*YL/age_per)*age_per

  dis_pri <- c()
  dis_pri_value <- mpl$dis_pri

  dis_pri[1:2] <- dis_pri_value
  dis1 <- dis_pri

  rho1 <- c()
  rho_prim_value <- mpl$rho_prim
  vacceff_prim <- mpl$vacceff_prim

  rho1[1] <- rho_prim_value
  rho1[2] <- rho_prim_value * (1 - vacceff_prim)

  phi_prim <- c()
  phi_prim_value <- mpl$phi_prim

  phi_prim[1:2] <- phi_prim_value
  phi1 <- phi_prim

  other_prop_immune <- mpl$other_prop_immune
  other_foi <- mpl$other_foi

  # initial compartment states

  aa <- bb <- c()

  for (i in 2:(na+1)){

    aa[i-1] <- (surv[i-1]-surv[i])/death[i-1]*(1-other_prop_immune*((1-exp(-other_foi*mean_age[i-1]))^1))

    bb[i-1] <- (surv[i-1]-surv[i])/death[i-1]*other_prop_immune*((1-exp(-other_foi*mean_age[i-1]))^1)

  }

  init_S <- array(0, c(na, 2, NP))

  init_S[1:na, 1, 1:NP] <- trunc(aa %o% Nb) # outer product

  init_I1 <- array(0, c(na, 2, NP))

  init_R1 <- array(0, c(na, 2, NP))

  init_R1[1:na, 1, 1:NP] <- trunc(bb %o% Nb)

  pTG <- mpl$propTransGlobal

  pTG_bigpatch <- pTG / 10

  res <- list(nn = nn,
              na = na,
              agec = agec,
              death = death,
              surv = surv,
              mean_surv = mean_surv,
              mean_age = mean_age,
              av_age_s = av_age_s,
              deathrt = deathrt,
              lifeexpec = lifeexpec,
              lifespan = lifespan,
              N_eq = N_eq,
              Nb = Nb,
              Mwt = init_Mwt_base,
              season_phase = season_phase,
              season_amp = season_amp,
              Wb_introtime = Wb_introtime,
              vacc_cu_rndtime = vacc_cu_rndtime,
              dis1 = dis1,
              rho1 = rho1,
              phi1 = phi1,
              init_S = init_S,
              init_I1 = init_I1,
              init_R1 = init_R1,
              pTG_bigpatch = pTG_bigpatch,
              pi = pi)

  append(res, mpl)

}
