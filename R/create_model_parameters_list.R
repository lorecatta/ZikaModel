#------------------------------------------------
#' model_param_list_create
#'
#' \code{model_param_list_create} creates a list of user-defined
#' parameters
#'
#' @export


model_param_list_create <- function(

                             # GENERAL

  YL = 364,                  # duration of year
  DT = 1,                    # size of time step
  NP = 21,                   # number of patches
  age_per = 52,              # time (weeks) between age updates
  N_human = 30000000,        # human population of each patch
  Omega = 1,
  DeltaBase = 0.2,           # adult mosquito mortality rate
  Sigma = 0.025,             # larval mosquito mortality rate
  Epsilon = 1/19,            # inverse of mean development time of larvae (1 / days)
  Rm = 2.69,                 # mosquito reproduction number (based on adult female fecundity and mortality rates)
  Mwt_mean = 1.5,            # mean of the lognorm distrib of eq number of adult female mosquitos per person, per patch, without seasonality
  Mwt_cv = 0.15,             # sd of the lognorm distrib ... (as above)
  eip_mean = 8.4,            # mean EIP
  incub = 5,                 # intrinsic incubation period (days)
  inf_per = 5,               # total duration of human infectiousness (days)
  nu = 1/21,                 # inverse of virus generation time (1 / days). Virus generation time = serial interval.
  Kappa = 0.5,               # biting rate per mosquito
  Beta_hm_1 = 0.7,           # per bite transmission probability from humans to mosquitoes
  Beta_mh_1 = 0.7,           # per bite transmission probability from mosquitoes to humans
  Delta_season = 1,          # scaling factor for effect of seasonality on adult mortality
  Kc_season = 1,             # scaling factor for effect of seasonality on larval carrying capacity
  eip_season = 1,            # scaling factor for effect of seasonality on EIP

                             # INTERVENTIONS

  propMwtControl = 0,        # Vector control / spraying (0.33)
  TimeMwtControlOn = 4,
  TimeMwtControlOff = 5,

  Wb_cyto = 1,               # Wolbachia
  Wb_mat = 1,
  Wb_fM = 0.95,
  Wb_fF = 0.95,
  Wb_relsusc1 = 0.9,
  Wb_relinf1 = 0.75,
  Wb_starttime = 140,
  Wb_introlevel = 0.5,
  Wb_introduration = 60,

  vacc_cu_minage = 2,        # Vaccine
  vacc_cu_maxage = 15,
  vacc_cu_coverage = 0.7,
  vacc_cu_time = 100,
  vacc_child_age = 3,
  vacc_child_coverage = 0.75,
  vacc_child_starttime = 30,
  vacc_child_stoptime = 30,

  dis_pri = 0.2,             # proportion of infections which are symptomatic
  rho_prim = 1,
  phi_prim = 1,
  vacceff_prim = 0.3,
  other_foi = 0.05,
  other_prop_immune = 0.5,
  propTransGlobal = 0.0005,  # proportion of transmission between all patches
  propTransNN = 0.005,       # proportion of transmssion with nearest-neighbor patches
  BG_FOI = 1e-6,
  AGE_REC = 2,
  PropDiseaseReported = 0.1

){

  # set up model param list
  mp_list <- list()

  ## DEFAULT PARAMETERS

  mp_list$YL <- YL
  mp_list$DT <- DT

  mp_list$NP <- NP
  mp_list$age_per <- age_per
  mp_list$N_human <- N_human
  mp_list$Omega <- Omega
  mp_list$DeltaBase <- DeltaBase
  mp_list$Sigma <- Sigma
  mp_list$Epsilon <- Epsilon
  mp_list$Rm <- Rm
  mp_list$Mwt_mean <- Mwt_mean
  mp_list$Mwt_cv <- Mwt_cv
  mp_list$eip_mean <- eip_mean
  mp_list$incub <- incub
  mp_list$inf_per <- inf_per
  mp_list$nu <- nu
  mp_list$Kappa <- Kappa
  mp_list$Beta_hm_1 <- Beta_hm_1
  mp_list$Beta_mh_1 <- Beta_mh_1
  mp_list$Delta_season <- Delta_season
  mp_list$Kc_season <- Kc_season
  mp_list$eip_season <- eip_season
  mp_list$propMwtControl <- propMwtControl
  mp_list$TimeMwtControlOn <- TimeMwtControlOn
  mp_list$TimeMwtControlOff <- TimeMwtControlOff
  mp_list$Wb_cyto <- Wb_cyto
  mp_list$Wb_mat <- Wb_mat
  mp_list$Wb_fM <- Wb_fM
  mp_list$Wb_fF <- Wb_fF
  mp_list$Wb_relsusc1 <- Wb_relsusc1
  mp_list$Wb_relinf1 <- Wb_relinf1
  mp_list$Wb_starttime <- Wb_starttime
  mp_list$Wb_introlevel <- Wb_introlevel
  mp_list$Wb_introduration <- Wb_introduration
  mp_list$vacc_cu_minage <- vacc_cu_minage
  mp_list$vacc_cu_maxage <- vacc_cu_maxage
  mp_list$vacc_cu_coverage <- vacc_cu_coverage
  mp_list$vacc_cu_time <- vacc_cu_time
  mp_list$vacc_child_age <- vacc_child_age
  mp_list$vacc_child_coverage <- vacc_child_coverage
  mp_list$vacc_child_starttime <- vacc_child_starttime
  mp_list$vacc_child_stoptime <- vacc_child_stoptime
  mp_list$dis_pri <- dis_pri
  mp_list$phi_prim <- phi_prim
  mp_list$rho_prim <- rho_prim
  mp_list$vacceff_prim <- vacceff_prim
  mp_list$other_foi <- other_foi
  mp_list$other_prop_immune <- other_prop_immune
  mp_list$propTransGlobal <- propTransGlobal
  mp_list$propTransNN <- propTransNN
  mp_list$BG_FOI <- BG_FOI
  mp_list$AGE_REC <- AGE_REC
  mp_list$PropDiseaseReported <- PropDiseaseReported

  mp_list

}
