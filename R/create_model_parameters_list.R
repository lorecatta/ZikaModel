#------------------------------------------------
#' model_param_list_create
#'
#' \code{model_param_list_create} creates a list of user-defined
#' parameters.
#'
#' @param YL Duration of a calendar year. Default = 364.
#' @param DT Time step size. Default = 1.
#' @param NP Number of patches. Default = 21.
#' @param age_per Time (weeks) between age updates. Default = 52 weeks.
#' @param N_human Human population in each patch. Default = 30000000.
#' @param Omega Density dependent power. Default = 1.
#' @param DeltaBase Adult mosquito mortality rate. Default = 0.2.
#' @param Sigma Larval mosquito mortality rate. Default = 0.025
#' @param Epsilon Larval mean development rate (1 / larvae mean development time in days). Default = 1/19.
#' @param Rm Mosquito reproduction number (based on adult female fecundity and mortality rates). Default = 2.69.
#' @param Mwt_mean Mean of the lognorm distribution of the at-equilibrium number of adult female mosquitos per person, per patch, without seasonality. Default = 1.5.
#' @param Mwt_cv Standard deviation of the lognorm distribution of the at-equilibrium number of adult female mosquitos per person, per patch, without seasonality. Default = 0.15
#' @param eip_mean Mean Extrinsic Incubation Period (days). Default = 8.4.
#' @param incub Intrinsic Incubation Period (days). Default = 5.5.
#' @param inf_per Total duration of human infectiousness (days). Default = 6.
#' @param nu Inverse of virus generation time (1 / days). Virus generation time = serial interval. Default = 1 / 21.
#' @param Kappa Mosquito biting rate. Default = 0.5.
#' @param Beta_hm_1 Per bite transmission probability from humans to mosquitoes. Default = 0.7.
#' @param Beta_mh_1 Per bite transmission probability from mosquitoes to humans. Default = 0.7.
#' @param season Logical for controlling the effect of seasonality.
#' @param propMwtControl Proportion of adult wild type mosquitoes which are killed by control. Default = 0.
#' @param TimeMwtControlOn Year of start of control of adult wild type mosquitoes. Default = 4.
#' @param TimeMwtControlOff Year of end of control of adult wild type mosquitoes. Default = 5.
#' @param Wb_cyto Default = 1.
#' @param Wb_mat Default = 1.
#' @param Wb_fM Default = 0.95.
#' @param Wb_fF Default = 0.95.
#' @param Wb_relsusc1 Default = 0.9.
#' @param Wb_relinf1 Default = 0.75.
#' @param Wb_starttime Default = 140.
#' @param Wb_introlevel Default = 0.5.
#' @param Wb_introduration Default = 60.
#' @param vacc_cu_minage Default = 2.
#' @param vacc_cu_maxage Default = 15.
#' @param vacc_cu_coverage Default = 0.7.
#' @param vacc_cu_time Default = 100.
#' @param vacc_child_age Deafult = 3.
#' @param vacc_child_coverage Default = 0.75.
#' @param vacc_child_starttime Default = 30.
#' @param vacc_child_stoptime Default = 30.
#' @param dis_pri Proportion of infections which are symptomatic. Default = 0.2.
#' @param rho_prim Default = 1.
#' @param phi_prim Default = 1.
#' @param vacceff_prim Efficacy of vaccination in reducing infection. Default = 0.3.
#' @param other_foi Default = 0.025.
#' @param other_prop_immune Default = 0.
#' @param propTransGlobal Proportion of transmission between all patches. Default = 0.0005.
#' @param propTransNN Proportion of transmssion with nearest-neighbor patches. Default = 0.005.
#' @param BG_FOI Background FOI. Default = 10e-8.
#' @param AGE_REC Default = 2.
#' @param PropDiseaseReported Reporting rate of symptomatic cases. Default = 0.1.

#' @export


model_param_list_create <- function(

  YL = 364,
  DT = 0.5,
  NP = 21,
  age_per = 1,
  N_human = 30000000,
  Omega = 1,
  DeltaBase = 0.2,
  Sigma = 0.025,
  Epsilon = 1/19,
  Rm = 2.69,
  Mwt_mean = 1.5,
  Mwt_cv = 0.15,
  eip_mean = 8.4,
  incub = 5.5,
  inf_per = 6,
  nu = 1/21,
  Kappa = 0.5,
  Beta_hm_1 = 0.7,
  Beta_mh_1 = 0.7,
  season,

  propMwtControl = 0,
  TimeMwtControlOn = 1.5,
  TimeMwtControlOff = 2.5,

  Wb_cyto = 1,
  Wb_mat = 1,
  Wb_fM = 0.95,
  Wb_fF = 0.95,
  Wb_relsusc1 = 0.9,
  Wb_relinf1 = 0.5,
  Wb_starttime = 140,
  Wb_introlevel = 0.5,
  Wb_introduration = 60,

  vacc_cu_minage = 2,
  vacc_cu_maxage = 15,
  vacc_cu_coverage = 0.7,
  vacc_cu_time = 150,
  vacc_child_age = 1,
  vacc_child_coverage = 0.56,
  vacc_child_starttime = 150,
  vacc_child_stoptime = 150,

  dis_pri = 0.2,
  rho_prim = 1,
  phi_prim = 1,
  vacceff_prim = 1,
  other_foi = 0.025,
  other_prop_immune = 0,
  propTransGlobal = 0.0005,
  propTransNN = 0.005,
  BG_FOI = 1e-8,
  AGE_REC = 2,
  PropDiseaseReported = 1){

  # Scaling factors (between 0 and 1) for effect of seasonality.
  # 1 = maximum effect of seasonality.
  # 0 = no effect of seasonality.

  if(season) {

    Kc_season <- 0.25
    eip_season <- 0.25
    Delta_season <- 0.25

  } else {

    Kc_season <- 0
    eip_season <- 0
    Delta_season <- 0

  }

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
