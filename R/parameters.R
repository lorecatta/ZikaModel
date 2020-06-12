
# -----------------------------------------------------------------------------

#' The function creates a list of user-defined parameters.
#'
#' @title Create a list of parameters
#'
#' @param YL Duration of a calendar year. Default = 365.
#' @param DT Time step size. Default = 1.
#' @param NP Number of patches. Default = 21.
#' @param agec Numeric vector of age group widths.
#'   Default = c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10).
#' @param death Numeric vector of mortality rates.
#'   Default = deathrt <- c(1e-10, 1e-10, 1e-10,
#'   0.00277068683332695, 0.0210680857689784, 0.026724997685722, 0.0525354529367476,
#'   0.0668013582441452, 0.119271483740379, 0.279105747097929, 0.390197266957464).
#' @param season Logical for controlling the effect of seasonality.
#'   TRUE = maximum effect of seasonality.
#'   FALSE = no effect of seasonality. Default = FALSE
#' @param age_per Time (weeks) between age updates. Default = 1.
#' @param N_human Human population in each patch. Default = 30000000.
#' #' @param incub Intrinsic Incubation Period (days). Default = 5.5.
#' @param inf_per Total duration of human infectiousness (days). Default = 6.
#' @param nu Inverse of virus generation time (1 / days).
#'   Virus generation time = serial interval. Default = 1 / 21.
#' @param Omega Intensity of density dependence of mosquito larvae mortality rate.
#'   0 = no density dependence. Default = 1.
#' @param DeltaBase Adult mosquito mortality rate. Default = 0.2.
#' @param Sigma Larval mosquito mortality rate. Default = 0.025
#' @param Epsilon Larval mean development rate (1 / larvae mean development time in days).
#'   Default = 1/19.
#' @param Rm Mosquito reproduction number (based on adult female fecundity and
#'   mortality rates). Default = 2.69.
#' @param Mwt_mean Mean of the lognorm distribution of the at-equilibrium number
#'   of adult female mosquitos per person, per patch, without seasonality.
#'   Default = 1.5.
#' @param Mwt_cv Standard deviation of the lognorm distribution of the
#'   at-equilibrium number of adult female mosquitos per person, per patch,
#'   without seasonality. Default = 0.15
#' @param eip_mean Mean Extrinsic Incubation Period (days). Default = 8.4.
#' @param Kappa Mosquito biting rate. Default = 0.5.
#' @param Beta_hm_1 Per bite transmission probability from humans to mosquitoes.
#'   Default = 0.7. This value is assigned to give a mean reproduction number, R0,
#'   across patches of 2.3 (with seasonal forcing).
#' @param Beta_mh_1 Per bite transmission probability from mosquitoes to humans.
#'   Arbitrarily assigned. Default = 0.7.
#' @param propMwtControl Increase in mortality of adult wild type mosquitoes induced
#'   by a \emph{general} type of intervention. Default = 0.
#' @param TimeMwtControlOn Year of start of control of adult wild type mosquitoes.
#'   Default = 2.
#' @param TimeMwtControlOff Year of end of control of adult wild type mosquitoes.
#'   Default = 3.
#' @param Wb_cyto Degree of cytoplasmic incompatibility induced by Wolbachia.
#'   Default = 1.
#' @param Wb_mat Degree of vertical transmission of Wolbachia. Default = 1.
#' @param Wb_fM Increase in mortality induced by Wolbachia. Default = 0.95.
#' @param Wb_fF Reduction in fecundity induced by Wolbachia. Default = 0.95.
#' @param Wb_relsusc1 Infectivity of a human host to Wolbachia infected mosquitoes
#'   (time τ after host infection). Default = 0.9.
#' @param Wb_relinf1 Infectiousness of Wolbachia-infected mosquitoes time τ after
#'   infection. Default = 0.75.
#' @param Wb_starttime Time of first release of Wolbachia-infected mosquiotes (years).
#'   Default = 1.
#' @param Wb_introlevel Ratio of Wolbachia-infected to wild type mosquitoes at
#'   introduction. (Not the proportion of Wolbachia AFTER introduction). Default = 0.
#' @param Wb_introduration Duration of Wolbachia release (days). Default = 60.
#' @param vacc_child_coverage Proportion of children vaccinated. Default = 0.
#' @param vacc_child_starttime Time when vaccination starts. Default = 30.
#' @param vacc_child_stoptime Time when vaccination stops. Default = 30.
#' @param vaccine_child_age Vector of binary indicators for routine vaccination
#'  of age groups. 1 = vaccinate, 0 = do not vaccinate. Same length as \code{agec}.
#'  Default = NULL.
#' @param vacc_cu_minage Minimum age at which children who missed vaccination can
#'   catch up. Default = 2.
#' @param vacc_cu_maxage Maximum age at which children who missed vaccination can
#'   catch up. Default = 15.
#' @param vacc_cu_coverage Proportion of children who undergo catch up vaccination.
#'   Default = 0.
#' @param vacc_cu_time Time when catch up vaccination occurs. Default = 30.
#' @param vaccine_cu_age Vector of binary indicators for catch up vaccination
#'  of age groups. 1 = vaccinate, 0 = do not vaccinate. Same length as \code{agec}.
#'  Default = c(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0), from 9 to 49.
#' @param vacceff_prim Efficacy of vaccination in reducing infection. Default = 0.75.
#' #' @param other_foi Default = 0.025.
#' @param other_prop_immune Propotion of population with pre-existing immunity.
#'   Default = 0.
#' @param propTransGlobal Proportion of transmission between all patches.
#'   Default = 0.0005.
#' @param propTransNN Proportion of transmssion with nearest-neighbor patches.
#'   Default = 0.
#' @param BG_FOI Force of infection on humans resulting from imported cases in
#'   travelers visiting from elsewhere. Default = 1e-8.
#' @param dis_pri Proportion of infections which are symptomatic. Default = 0.2.
#' @param rho_prim Default = 1.
#' @param phi_prim Default = 1.
#' @param AGE_REC Default = 2.
#' @param PropDiseaseReported Reporting rate of symptomatic cases. Default = 0.1.
#'
#' @return Parameter list
#'
#' @export
parameters_deterministic_model <- function(

  # initial state, duration, patches
  YL = 365,
  DT = 1,
  NP = 21,

  # demography
  agec = default_demog$agec,
  death = default_demog$death,

  # seasonality
  season = FALSE,

  # humans
  age_per = 1,
  N_human = 30000000,
  incub = 5.5,
  inf_per = 6,
  nu = 1/21,

  # mosquitoes
  Omega = 1,
  DeltaBase = 0.2,
  Sigma = 0.025,
  Epsilon = 1/19,
  Rm = 2.69,
  Mwt_mean = 1.5,
  Mwt_cv = 0.15,
  eip_mean = 8.4,
  Kappa = 0.5,
  Beta_hm_1 = 0.7,
  Beta_mh_1 = 0.7,

  # general control
  propMwtControl = 0,
  TimeMwtControlOn = 1.5,
  TimeMwtControlOff = 2.5,

  # wolbachia
  Wb_cyto = 1,
  Wb_mat = 1,
  Wb_fM = 0.95,
  Wb_fF = 0.95,
  Wb_relsusc1 = 0.9,
  Wb_relinf1 = 0.5,
  Wb_starttime = 1,
  Wb_introlevel = 0,
  Wb_introduration = 60,

  # vaccination
  vacc_child_coverage = 0,
  vacc_child_starttime = 30,
  vacc_child_stoptime = 30,
  vaccine_child_age = NULL,
  vacc_cu_minage = 2,
  vacc_cu_maxage = 15,
  vacc_cu_coverage = 0,
  vacc_cu_time = 30,
  vaccine_cu_age = c(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
  vacceff_prim = 0.75,

  # FOI
  other_foi = 0.025,
  other_prop_immune = 0,
  propTransGlobal = 0.0005,
  propTransNN = 0,
  BG_FOI = 1e-8,

  # disease
  dis_pri = 0.2,
  rho_prim = 1,
  phi_prim = 1,
  AGE_REC = 2,
  PropDiseaseReported = 1) {

  # Check parameters
  if(!is.null(vaccine_child_age)) {
    assert_same_length(vaccine_child_age, agec)
  }
  assert_same_length(vaccine_cu_age, agec)

  pars_to_equlibrium_init_create <- list(NP = NP,
                                         YL = YL,
                                         DT = DT,
                                         N_human = N_human,
                                         other_prop_immune = other_prop_immune,
                                         other_foi = other_foi)

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec,
                                        death,
                                        pars_to_equlibrium_init_create)

  amp_phas <- ZikaModel::amplitudes_phases # Amplitude and phase of seasonal forcing for each patch.
  amp_phas[, "phase"] <- amp_phas[, "phase"] * YL

  Wb_introtime <- c()

  for (i in seq_len(NP - 1)){

    Wb_introtime[i] <- Wb_starttime # + trunc(6 * (i - 1) / NP) * 0.5

  }

  Wb_introtime[NP] <- 1000

  vacc_cu_rndtime <- trunc(vacc_cu_time*YL/age_per)*age_per

  dis_pri_values <- c()
  dis_pri_values[1:2] <- dis_pri

  rho1_values <- c()
  rho1_values[1] <- rho_prim
  rho1_values[2] <- rho_prim * (1 - vacceff_prim)

  phi_prim_values <- c()
  phi_prim_values[1:2] <- phi_prim

  pTG_bigpatch <- propTransGlobal / 10

  na <- as.integer(length(agec))  # number of age groups

  if(!is.null(vaccine_child_age)) {

    vaccine_child_age_2 <- c(0, vaccine_child_age)

  } else {

    vaccine_child_age_2 <- rep(0, na + 1)

  }

  if(!is.null(vaccine_cu_age)) {

    vaccine_cu_age_2 <- c(0, vaccine_cu_age)

  } else {

    vaccine_cu_age_2 <- rep(0, na + 1)

  }

  # Scaling factors (between 0 and 1) for effect of seasonality.
  # 1 = maximum effect of seasonality.
  # 0 = no effect of seasonality.

  if(season) {

    Kc_season <- 1
    eip_season <- 1
    Delta_season <- 1

  } else {

    Kc_season <- 0
    eip_season <- 0
    Delta_season <- 0

  }

  # Collate Parameters Into List
  mp_list <- list(YL = YL,
                  DT = DT,
                  NP = NP,
                  agec = agec,
                  death = death,
                  age_per = age_per,
                  N_human = N_human,
                  incub = incub,
                  inf_per = inf_per,
                  nu = nu,
                  Omega = Omega,
                  DeltaBase = DeltaBase,
                  Sigma = Sigma,
                  Epsilon = Epsilon,
                  Rm = Rm,
                  Mwt_mean = Mwt_mean,
                  Mwt_cv = Mwt_cv,
                  eip_mean = eip_mean,
                  Kappa = Kappa,
                  Beta_hm_1 = Beta_hm_1,
                  Beta_mh_1 = Beta_mh_1,
                  propMwtControl = propMwtControl,
                  TimeMwtControlOn = TimeMwtControlOn,
                  TimeMwtControlOff = TimeMwtControlOff,
                  Wb_cyto = Wb_cyto,
                  Wb_mat = Wb_mat,
                  Wb_fM = Wb_fM,
                  Wb_fF = Wb_fF,
                  Wb_relsusc1 = Wb_relsusc1,
                  Wb_relinf1 = Wb_relinf1,
                  Wb_starttime = Wb_starttime,
                  Wb_introlevel = Wb_introlevel,
                  Wb_introduration = Wb_introduration,
                  vacc_child_coverage = vacc_child_coverage,
                  vacc_child_starttime = vacc_child_starttime,
                  vacc_child_stoptime = vacc_child_stoptime,
                  vacc_cu_minage = vacc_cu_minage,
                  vacc_cu_maxage = vacc_cu_maxage,
                  vacc_cu_coverage = vacc_cu_coverage,
                  vacc_cu_time = vacc_cu_time,
                  propTransGlobal = propTransGlobal,
                  propTransNN = propTransNN,
                  BG_FOI = BG_FOI,
                  dis1 = dis_pri_values,
                  rho1 = rho1_values,
                  phi1 = phi_prim_values,
                  AGE_REC = AGE_REC,
                  PropDiseaseReported = PropDiseaseReported,
                  Wb_introtime = Wb_introtime,
                  vacc_cu_rndtime = vacc_cu_rndtime,
                  nn = state_init$nn,
                  na = state_init$na,
                  surv = state_init$surv,
                  mean_surv = state_init$mean_surv,
                  mean_age = state_init$mean_age,
                  av_age_s = state_init$av_age_s,
                  deathrt = state_init$deathrt,
                  lifeexpec = state_init$lifeexpec,
                  lifespan = state_init$lifespan,
                  N_eq = state_init$N_eq,
                  Nb = state_init$Nb,
                  init_S = state_init$init_S_sign,
                  init_I1 = state_init$init_I1,
                  init_R1 = state_init$init_R1,
                  pTG_bigpatch = pTG_bigpatch,
                  vacc_child_age = vaccine_child_age_2,
                  vacc_cu_age = vaccine_cu_age_2,
                  Kc_season = Kc_season,
                  eip_season = eip_season,
                  Delta_season = Delta_season,
                  pi = pi)

  return(mp_list)

}
