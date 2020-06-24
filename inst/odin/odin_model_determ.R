
# -----------------------------------------------------------------------------
#
# setting parameters for the model
#
# -----------------------------------------------------------------------------


# general ---------------------------------------------------------------------


DT <- user() # time step
time <- step * DT # time value
YL <- user() # year length
pi <- user()

# number of patches
NP <- user()
# 8 nearest-neigbors to each patch
nn[,] <- user()

# amplitude and phases of the seasonal forcing
amplitudes_phases[,] <- user()
# scaling factor for effect of seasonality on adult mosquito mortality rate
Delta_season <- user()
# scaling factor for effect of seasonality on mosquito larvae carrying capacity
Kc_season <- user()
# scaling factor for effect of seasonality on extrinsic incubation period
eip_season <- user()


# mosquitoes ------------------------------------------------------------------


# mosquito reproduction number
Rm <- user()
# larvae development rate
Epsilon <- user()
# larvae mortality rate
Sigma <- user()
# Intensity of density dependence of mosquito larvae mortality rate
Omega <- user()
# biting rate
Kappa <- user()
# mean of the lognorm distribution of number of adult female mosquitos per person,
# per patch, without seasonality, at equlibrium.
Mwt_mean <- user()
# standard deviation of the lognorm distribution
Mwt_cv <- user()
# adult mosquito mortality rate
DeltaBase <- user()
# mean extrinsic incubation period
eip_mean <- user()
# per bite transmission probability from mosquitoes to humans
Beta_mh_1 <- user()
# per bite transmission probability from humans to mosquitoes
Beta_hm_1 <- user()


# humans ----------------------------------------------------------------------


# intrinsic incubation period
incub <- user()
# duration of human infectiousness
inf_per <- user()
# recovery rate
nu <- user()
# number of births
Nb[] <- user()
# human population in each patch
N_eq[] <- user()
# mortality rates
deathrt[] <- user()
# time (weeks) between age updates
age_per <- user()
# number of age groups
na <- user()
# width of each age group
agec[] <- user()


# interventions ---------------------------------------------------------------


# year of start of control of adult wild type mosquitoes
TimeMwtControlOn <- user()
# year of end of control of adult wild type mosquitoes
TimeMwtControlOff <- user()
# increase in mortality of adult wild type mosquitoes induced by general control
propMwtControl <- user()

# cytoplasmic incompatibility induced by Wolbachia
Wb_cyto <- user()
# degree of vertical transmission of Wolbachia
Wb_mat <- user()
# increase in mortality induced by Wolbachia
Wb_fM <- user()
# reduction in fecundity induced by Wolbachia
Wb_fF <- user()

# Infectivity of a human host to Wolbachia infected mosquitoes (time τ after host infection).
Wb_relsusc1 <- user()
# infectiousness of Wolbachia-infected mosquitoes time τ after infection
Wb_relinf1 <- user()
# time of first release of Wolbachia-infected mosquiotes (years)
Wb_introtime[] <- user()
# ratio of Wolbachia-infected to wild type mosquitoes at introduction
Wb_introlevel <- user()
# duration of Wolbachia release (days)
Wb_introduration <- user()

# decay of protection from vaccine
rho1[] <- user()
# proportion of age group undergoing routine vaccination
vacc_child_coverage <- user()
# time when routine vaccination starts
vacc_child_starttime <- user()
# time when routine vaccination ends
vacc_child_stoptime <- user()
# vector of binary indicators indicating which age groups undergo routine vaccination
vacc_child_age[] <- user()

# proportion of age group undergoing catch up vaccination
vacc_cu_coverage <- user()
# time when catch up vaccination occurs
vacc_cu_rndtime <- user()
# vector of binary indicators indicating which age groups undergo catch up vaccination
vacc_cu_age[] <- user()


# FOI -------------------------------------------------------------------------


# proportion of transmission between all patches.
propTransGlobal <- user()
# proportion of transmission with nearest-neighbor patches
propTransNN <- user()
# force of infection on humans resulting from imported cases in travelers
# visiting from elsewhere
BG_FOI <- user()
# degree of cross-immunity  (not applicable for Zika)
phi1[] <- user()


# initial arrays --------------------------------------------------------------


# susceptibles
init_S[,,] <- user()
# exposed and infectious
init_I1[,,] <- user()
# recovered
init_R1[,,] <- user()


# -----------------------------------------------------------------------------


vnc_row <- na + 1 # position 1 is like position 0 in BM

MwtCont <- if ((time >= TimeMwtControlOn * YL) &&
               (time < TimeMwtControlOff * YL)) (1 - propMwtControl) else 1

# Number of eggs produced per female mosquito per time step (baseline fecundity). Rm is fixed.
Gamma <- Rm * DeltaBase * (Epsilon + Sigma) / Epsilon


ln_sd <- sqrt(log(1 + Mwt_cv * Mwt_cv))
ln_mean <- log(Mwt_mean) - (ln_sd * ln_sd) / 2

dim(ln_pick) <- NP
ln_pick[] <- ln_mean

Mwt[] <- exp(ln_pick[i])


# -----------------------------------------------------------------------------
# Mean value of mosquito mortality, carrying capacity and EIP
# without the effect of seasonality


# mean adult mosquito mortality rate, corrected by the effect of control
DeltaMean <- DeltaBase / MwtCont

# mean larval mosquito carrying capacity. Mwt (in Kc equation below) is fixed.
Kc_mean <- DeltaBase *
  ((Epsilon * (Gamma - DeltaBase) / (DeltaBase * Sigma) - 1) ^ (-1 / Omega)) / Epsilon


# -----------------------------------------------------------------------------
# Mean value of mosquito mortality, carrying capacity and EIP
# with the effect of seasonality


season_amp[] <- amplitudes_phases[i, 2]
season_phases[] <- amplitudes_phases[i, 3]

Delta[1:(NP-1)] <- DeltaMean /
  (1 + season_amp[i] * Delta_season * cos(2 * pi * (time + season_phases[i]) / YL))
Delta[NP] <- DeltaMean

Kc[1:(NP-1)] <- Mwt[i] * Kc_mean *
  (1 + season_amp[i] * Kc_season * cos(2 * pi * (time + season_phases[i]) / YL))
Kc[NP] <- MwtCont * Mwt_mean * Kc_mean *
  (1 + season_amp[i] * Kc_season * cos(2 * pi * (time + season_phases[i]) / YL))

eip[] <- eip_mean *
  (1 - season_amp[i] * eip_season * cos(2 * pi * (time + season_phases[i]) / YL))

Delta_wb[] <- Delta[i] / Wb_fM

vacc_noncov[, 1] <- (if ((time >= YL * vacc_child_starttime) && (time < YL * vacc_child_stoptime) && (vacc_child_age[i] == 1)) (1 - vacc_child_coverage) else 1)
vacc_noncov[, 2] <- 1

vcu_noncov[, 1] <- (if ((time == vacc_cu_rndtime) && (vacc_cu_age[i] == 1)) (1 - vacc_cu_coverage) else 1)
vcu_noncov[, 2] <- 1



# -----------------------------------------------------------------------------
#
# States of wild type mosquitoes
#
# -----------------------------------------------------------------------------



Mwt_tot[] <- Mwt_S[i] + Mwt_E1[i] + Mwt_E2[i] + Mwt_I1[i]

# Fecundity (i.e., rate at which female mosquitoes lay eggs)
# with effect of cytoplasmic incompatibility
Lwt_birth_lambda[] <- DT *
  (Gamma * Mwt_tot[i] * (Mwt_tot[i] + (1 - Wb_cyto) * Mwb_tot[i]) /
     (Mwt_tot[i] + Mwb_tot[i]) + Wb_fF * (1 - Wb_mat) * Mwb_tot[i])
Lwt_birth[] <- Lwt_birth_lambda[i]

L_deathrt[] <- DT * Sigma * ((1 + ((Lwt[i] + Lwb[i]) / (Kc[i] * NTp[i])) ^ Omega))

L_dr[] <- if (L_deathrt[i] >= 1) 0.999 else L_deathrt[i]

O_Lwt_prob[] <- max(min(DT * Epsilon + L_dr[i], 1), 0)
O_Lwt[] <- Lwt[i] * O_Lwt_prob[i]

Lwt_mature_prob[] <- max(min(DT * Epsilon / (DT * Epsilon + L_dr[i]), 1), 0)
Lwt_mature[] <- O_Lwt[i] * Lwt_mature_prob[i]

update(Lwt[]) <- Lwt_birth[i] + Lwt[i] - O_Lwt[i]

Mwt_FOI1[] <- DT * Beta_hm_1 * Kappa * infectious1[i]

O_Mwt_S_prob[] <- max(min(DT * Delta[i] + Mwt_FOI1[i], 1), 0)
O_Mwt_S[] <- Mwt_S[i] * O_Mwt_S_prob[i]

Mwt_inf1_prob[] <- max(min(Mwt_FOI1[i] / (DT * Delta[i] + Mwt_FOI1[i]), 1), 0)
Mwt_inf1[] <- O_Mwt_S[i] * Mwt_inf1_prob[i]

update(Mwt_S[]) <- Lwt_mature[i] + Mwt_S[i] - O_Mwt_S[i]

O_Mwt_E1_prob[] <- max(min(DT * (Delta[i] + 2 / eip[i]), 1), 0)
O_Mwt_E1[] <- Mwt_E1[i] * O_Mwt_E1_prob[i]

Mwt_E1_incub_prob[] <- max(min(1 / (Delta[i] * eip[i] / 2 + 1), 1), 0)
Mwt_E1_incub[] <- O_Mwt_E1[i] * Mwt_E1_incub_prob[i]

update(Mwt_E1[]) <- Mwt_inf1[i] + Mwt_E1[i] - O_Mwt_E1[i]

O_Mwt_E2_prob[] <- max(min(DT * (Delta[i] + 2 / eip[i]), 1), 0)
O_Mwt_E2[] <- Mwt_E2[i] * O_Mwt_E2_prob[i]

Mwt_E2_incub_prob[] <- max(min(1 / (Delta[i] * eip[i] / 2 + 1), 1), 0)
Mwt_E2_incub[] <- O_Mwt_E2[i] * Mwt_E2_incub_prob[i]

update(Mwt_E2[]) <- Mwt_E1_incub[i] + Mwt_E2[i] - O_Mwt_E2[i]

O_Mwt_I1_prob[] <- max(min(DT * Delta[i], 1), 0)
O_Mwt_I1[] <- Mwt_I1[i] * O_Mwt_I1_prob[i]

update(Mwt_I1[]) <- Mwt_E2_incub[i] + Mwt_I1[i] - O_Mwt_I1[i]



# -----------------------------------------------------------------------------
#
# States of wolbachia type mosquitoes
#
# -----------------------------------------------------------------------------



Mwb_tot[] <- Mwb_S[i] + Mwb_E1[i] + Mwb_E2[i] + Mwb_I1[i]

R0t_1[] <- Kappa * Kappa * (Mwt_tot[i] + Wb_relsusc1 * Wb_relinf1 * Mwb_tot[i]) *
  Beta_hm_1 * inf_per * Beta_mh_1 / (1 + Delta[i] * eip[i]) / Delta[i] / NTp[i]

Lwb_birth_lambda[] <- DT *(Gamma * Wb_fF * Wb_mat * Mwb_tot[i])
Lwb_birth[] <- Lwb_birth_lambda[i]

O_Lwb_prob[] <- max(min(DT * Epsilon + L_dr[i], 1), 0)
O_Lwb[] <- Lwb[i] * O_Lwb_prob[i]

Lwb_mature_prob[] <- max(min(DT * Epsilon / (DT * Epsilon + L_dr[i]), 1), 0)
Lwb_mature[] <- O_Lwb[i] * Lwb_mature_prob[i]

update(Lwb[]) <- Lwb_birth[i] + Lwb[i] - O_Lwb[i]

Mwb_FOI1[] <- Wb_relsusc1 * Mwt_FOI1[i]

O_Mwb_S_prob[] <- max(min(DT * Delta_wb[i] + Mwb_FOI1[i], 1), 0)
O_Mwb_S[] <- Mwb_S[i] * O_Mwb_S_prob[i]

Mwb_inf1_prob[] <- max(min(Mwb_FOI1[i] / (DT * Delta_wb[i] + Mwb_FOI1[i]), 1), 0)
Mwb_inf1[] <- O_Mwb_S[i] * Mwb_inf1_prob[i]

Wb_introrate[] <- Wb_introlevel * Mwt[i] * N_eq[i] * DT / Wb_introduration

Mwb_intro[] <- if((time >= Wb_introtime[i] * YL) &&
                  (time < Wb_introtime[i] * YL + Wb_introduration)) Wb_introrate[i] else 0

update(Mwb_S[]) <- Lwb_mature[i] + Mwb_intro[i] + Mwb_S[i] - O_Mwb_S[i]

O_Mwb_E1_prob[] <- max(min(DT * (Delta_wb[i] + 1 / eip[i]), 1), 0)
O_Mwb_E1[] <- Mwb_E1[i] * O_Mwb_E1_prob[i]

Mwb_E1_incub_prob[] <- max(min(1 / (Delta[i] * eip[i] + 1), 1), 0)
Mwb_E1_incub[] <- O_Mwb_E1[i] * Mwb_E1_incub_prob[i]

update(Mwb_E1[]) <- Mwb_inf1[i] + Mwb_E1[i] - O_Mwb_E1[i]

O_Mwb_E2_prob[] <- max(min(DT * (Delta[i] + 2 / eip[i]), 1), 0)
O_Mwb_E2[] <- Mwb_E2[i] * O_Mwb_E2_prob[i]

Mwb_E2_incub_prob[] <- max(min(1 / (Delta[i] * eip[i] / 2 + 1), 1), 0)
Mwb_E2_incub[] <- O_Mwb_E2[i] * Mwb_E2_incub_prob[i]

update(Mwb_E2[]) <- Mwb_E1_incub[i] + Mwb_E2[i] - O_Mwb_E2[i]

O_Mwb_I1_prob[] <- max(min(DT * Delta_wb[i], 1), 0)
O_Mwb_I1[] <- Mwb_I1[i] * O_Mwb_I1_prob[i]

update(Mwb_I1[]) <- Mwb_E2_incub[i] + Mwb_I1[i] - O_Mwb_I1[i]



# -----------------------------------------------------------------------------
#
# Human states
#
# -----------------------------------------------------------------------------



agerts <- if (trunc(time / age_per) == time / age_per) age_per / YL else 0
# agerts=DT/YL
agert[] <- agerts / agec[i] # ageing rate per age group


births_lambda[] <- DT * Nb[i] / YL
births[] <- births_lambda[i]

O_S_prob[,,] <- max(min(rho1[j] * FOI1[k] + agert[i] + deathrt[i], 1), 0)
O_S[,,] <- S[i,j,k] * O_S_prob[i,j,k]

inf_1_prob[,,] <- max(min(rho1[j] * FOI1[k] / (rho1[j] * FOI1[k] +
                                                 agert[i] + deathrt[i]), 1), 0)
inf_1[,,] <- O_S[i,j,k] * inf_1_prob[i,j,k]

age_S_trials[,,] <- O_S[i,j,k] - inf_1[i,j,k]
age_S_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_S[,,] <- age_S_trials[i,j,k] * age_S_prob[i]

dim(n_S) <- c(na, 2, NP)
n_S[1,1,1:NP] <- trunc(0.5 + births[k] + S[i,j,k] - O_S[i,j,k])
n_S[2:na,1,1:NP] <- trunc(0.5 + vacc_noncov[i,j] * age_S[i-1,j,k] + (1 - vacc_noncov[i,3-j]) * age_S[i-1,3-j,k] + S[i,j,k] - O_S[i,j,k])

n_S[1,2,1:NP] <- trunc(0.5 + S[i,j,k] - O_S[i,j,k])
n_S[2:na,2,1:NP] <- trunc(0.5 + vacc_noncov[i,j] * age_S[i-1,j,k] + (1 - vacc_noncov[i,3-j]) * age_S[i-1,3-j,k] + S[i,j,k] - O_S[i,j,k])

update(S[1:na,1,1:NP]) <- vcu_noncov[i,j] * n_S[i,j,k]
update(S[1:na,2,1:NP]) <- (1 - vcu_noncov[i,1]) * n_S[i,1,k] + n_S[i,2,k]

O_I1_prob[] <- max(min(nu + agert[i] + deathrt[i], 1), 0)
O_I1[,,] <- I1[i,j,k] * O_I1_prob[i]

recov1_prob[] <- max(min(nu / (nu + agert[i] + deathrt[i]), 1), 0)
recov1[,,] <- O_I1[i,j,k] * recov1_prob[i]

age_I1_trials[,,] <- O_I1[i,j,k] - recov1[i,j,k]
age_I1_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_I1[2:vnc_row,1:2,1:NP] <- age_I1_trials[i-1,j,k] * age_I1_prob[i-1]
age_I1[1,1:2,1:NP] <- 0

dim(n_I1) <- c(na, 2, NP)
n_I1[1:na,1:2,1:NP] <- trunc(0.5 + vacc_noncov[i,j] * age_I1[i,j,k] + (1 - vacc_noncov[i,3-j]) * age_I1[i,3-j,k] + inf_1[i,j,k] + I1[i,j,k] - O_I1[i,j,k])

update(I1[1:na,1,1:NP]) <- vcu_noncov[i,j] * n_I1[i,j,k]
update(I1[1:na,2,1:NP]) <- (1 - vcu_noncov[i,1]) * n_I1[i,1,k] + n_I1[i,2,k]

O_R1_prob[] <- max(min(agert[i] + deathrt[i], 1), 0)
O_R1[,,] <- R1[i,j,k] * O_R1_prob[i]

age_R1_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_R1[2:vnc_row,1:2,1:NP] <- O_R1[i-1,j,k] * age_R1_prob[i-1]
age_R1[1,1:2,1:NP] <- 0

dim(n_R1) <- c(na, 2, NP)
n_R1[1:na,1:2,1:NP] <- trunc(0.5 + vacc_noncov[i,j] * age_R1[i,j,k] + (1 - vacc_noncov[i,3-j]) * age_R1[i,3-j,k] + recov1[i,j,k] + R1[i,j,k] - O_R1[i,j,k])

update(R1[1:na,1,1:NP]) <- vcu_noncov[i,j] * n_R1[i,j,k]
update(R1[1:na,2,1:NP]) <- (1 - vcu_noncov[i,1]) * n_R1[i,1,k] + n_R1[i,2,k]



# -----------------------------------------------------------------------------
#
# Keep track of infectious humans
#
# -----------------------------------------------------------------------------



# 2 human incubation and infectious classes (instead of 1)
# are needed for a gamma rather than
# an exponentially distributed human infectious period


Ntotal[,,] <- S[i,j,k] + I1[i,j,k] + R1[i,j,k]

Ntotal_sum[1,1:NP] <- Ntotal[i,1,j] + Ntotal[i,2,j]
Ntotal_sum[2:na,1:NP] <- Ntotal_sum[i-1,j] + Ntotal[i,1,j] + Ntotal[i,2,j]
NTp[] <- Ntotal_sum[na,i]

Y1[,,] <- phi1[j] * inf_1[i,j,k]

Y1T_sum[1,1:NP] <- Y1[i,1,j] + Y1[i,2,j]
Y1T_sum[2:na,1:NP] <- Y1T_sum[i-1,j] + Y1[i,1,j] + Y1[i,2,j]

Y1T[] <- Y1T_sum[na,i] / NTp[i]

update(incubA[]) <- incubA[i] + Y1T[i] - DT * 2 * incubA[i] / incub

update(incubB[]) <- incubB[i] + DT * 2 * (incubA[i] - incubB[i]) / incub

update(infectiousA[]) <- infectiousA[i] +
  DT * 2 * (incubB[i] / incub - infectiousA[i] / inf_per)

update(infectiousB[]) <- infectiousB[i] +
  DT * 2 * (infectiousA[i] - infectiousB[i]) / inf_per

infectious1[] <- infectiousA[i] + infectiousB[i]

# initial(infectious1[]) <- 0
# update(infectious1[]) <- Y1T_del_inc[i] + infectious1[i] - Y1T_del_inc_ip[i]
# incub_2 <- incub / DT
# Y1T_del_inc[] <- delay(Y1T[i], incub_2)
# inf_per_2 <- inf_per / DT
# Y1T_lag <- incub_2 + inf_per_2
# Y1T_del_inc_ip[] <- delay(Y1T[i], Y1T_lag)



# -----------------------------------------------------------------------------
#
# Calculate FOI
#
# -----------------------------------------------------------------------------



FOI1p[] <- DT * Beta_mh_1 * Kappa * (Mwt_I1[i] + Mwb_I1[i] * Wb_relinf1) / NTp[i]

FOI1nn[1:(NP-1)] <- (FOI1p[nn[i,1]] + FOI1p[nn[i,2]] + FOI1p[nn[i,3]] +
                     FOI1p[nn[i,4]] + FOI1p[nn[i,5]] + FOI1p[nn[i,6]] +
                     FOI1p[nn[i,7]] + FOI1p[nn[i,8]]) / 8

FOI1av <- (sum(FOI1p[]) - FOI1p[NP]) / (NP - 1)

FOI1[1:(NP-1)] <- propTransNN * FOI1nn[i] + propTransGlobal * FOI1av +
  (1 - propTransGlobal - propTransNN) * FOI1p[i] + DT * BG_FOI / YL
# FOI1[NP] <- pTG_bigpatch * FOI1av + (1 - pTG_bigpatch) * FOI1p[i]
FOI1[NP] <- 0



# -----------------------------------------------------------------------------



## initial states
initial(Lwt[]) <- trunc(Mwt[i] * N_eq[i] * Delta[i] / Epsilon)
initial(Mwt_S[]) <- trunc(Mwt[i] * N_eq[i])
initial(Mwt_E1[]) <- 0
initial(Mwt_E2[]) <- 0
initial(Mwt_I1[]) <- 0
initial(Lwb[]) <- 0
initial(Mwb_S[]) <- 0
initial(Mwb_E1[]) <- 0
initial(Mwb_E2[]) <- 0
initial(Mwb_I1[]) <- 0
initial(S[,,]) <- init_S[i,j,k]
initial(I1[,,]) <- init_I1[i,j,k]
initial(R1[,,]) <- init_R1[i,j,k]
initial(incubA[]) <- 0
initial(incubB[]) <- 0
initial(infectiousA[]) <- 0
initial(infectiousB[]) <- 0

## array dimensions
# for the initial values
dim(init_S) <- c(na, 2, NP)
dim(init_I1) <- c(na, 2, NP)
dim(init_R1) <- c(na, 2, NP)

# for state variables
dim(Lwt) <- NP
dim(Mwt_S) <- NP
dim(Mwt_E1) <- NP
dim(Mwt_E2) <- NP
dim(Mwt_I1) <- NP
dim(Lwb) <- NP
dim(Mwb_S) <- NP
dim(Mwb_E1) <- NP
dim(Mwb_E2) <- NP
dim(Mwb_I1) <- NP
dim(S) <- c(na, 2, NP)
dim(I1) <- c(na, 2, NP)
dim(R1) <- c(na, 2, NP)

# for the number of mosquitoes / people moving In and Out of compartments
dim(Lwt_birth) <- NP
dim(O_Lwt) <- NP
dim(Lwt_mature) <- NP
dim(O_Mwt_S) <- NP
dim(Mwt_inf1) <- NP
dim(O_Mwt_E1) <- NP
dim(Mwt_E1_incub) <- NP
dim(O_Mwt_E2) <- NP
dim(Mwt_E2_incub) <- NP
dim(O_Mwt_I1) <- NP
dim(Lwb_birth) <- NP
dim(O_Lwb) <- NP
dim(Lwb_mature) <- NP
dim(O_Mwb_S) <- NP
dim(Mwb_inf1) <- NP
dim(O_Mwb_E1) <- NP
dim(Mwb_E1_incub) <- NP
dim(Mwb_E2_incub) <- NP
dim(O_Mwb_I1) <- NP
dim(births) <- NP
dim(O_S) <- c(na, 2, NP)
dim(inf_1) <- c(na, 2, NP)
dim(age_S) <- c(na, 2, NP)
dim(age_S_trials) <- c(na, 2, NP)
dim(O_I1) <- c(na, 2, NP)
dim(recov1) <- c(na, 2, NP)
dim(age_I1) <- c(vnc_row, 2, NP)
dim(age_I1_trials) <- c(na, 2, NP)
dim(O_R1) <- c(na, 2, NP)
dim(age_R1) <- c(vnc_row, 2, NP)
dim(incubA) <- NP
dim(incubB) <- NP
dim(infectiousA) <- NP
dim(infectiousB) <- NP
dim(infectious1) <- NP

# for the probabilities of generating the number of mosquitoes / people
# moving In and Out of compartments
dim(O_Lwt_prob) <- NP
dim(Lwt_mature_prob) <- NP
dim(O_Mwt_S_prob) <- NP
dim(Mwt_inf1_prob) <- NP
dim(O_Mwt_E1_prob) <- NP
dim(Mwt_E1_incub_prob) <- NP
dim(O_Mwt_E2_prob) <- NP
dim(Mwt_E2_incub_prob) <- NP
dim(O_Mwt_I1_prob) <- NP
dim(O_Lwb_prob) <- NP
dim(Lwb_mature_prob) <- NP
dim(O_Mwb_S_prob) <- NP
dim(Mwb_inf1_prob) <- NP
dim(O_Mwb_E1_prob) <- NP
dim(Mwb_E1_incub_prob) <- NP
dim(O_Mwb_E2_prob) <- NP
dim(O_Mwb_E2) <- NP
dim(Mwb_E2_incub_prob) <- NP
dim(O_Mwb_I1_prob) <- NP
dim(O_S_prob) <- c(na, 2, NP)
dim(inf_1_prob) <- c(na, 2, NP)
dim(age_S_prob) <- na
dim(O_I1_prob) <- na
dim(recov1_prob) <- na
dim(age_I1_prob) <- na
dim(O_R1_prob) <- na
dim(age_R1_prob) <- na

# for the rates
dim(Lwt_birth_lambda) <- NP
dim(Lwb_birth_lambda) <- NP
dim(L_deathrt) <- NP
dim(L_dr) <- NP
dim(Nb) <- NP
dim(births_lambda) <- NP
dim(deathrt) <- na
dim(agert) <- na

# for FOI/R0 calculation
dim(Mwt_FOI1) <- NP
dim(Mwb_FOI1) <- NP
dim(FOI1p) <- NP
dim(FOI1nn) <- NP
dim(FOI1) <- NP
dim(R0t_1) <- NP
dim(phi1) <- 2

# for intermediate outputs
dim(agec) <- na
dim(Mwt) <- NP
dim(Mwt_tot) <- NP
dim(Mwb_tot) <- NP
dim(N_eq) <- NP
dim(nn) <- c(NP, 8)
dim(amplitudes_phases) <- c(NP, 3)
dim(season_phases) <- NP
dim(season_amp) <- NP
dim(Delta) <- NP
dim(Delta_wb) <- NP
dim(Kc) <- NP
dim(eip) <- NP
dim(Ntotal) <- c(na, 2, NP)
dim(Ntotal_sum) <- c(na, NP)
dim(NTp) <- NP
dim(Y1) <- c(na, 2, NP)
dim(Y1T_sum) <- c(na,NP)
dim(Y1T) <- NP

# for the interventions
dim(Wb_introrate) <- NP
dim(Mwb_intro) <- NP
dim(Wb_introtime) <- NP
dim(rho1) <- 2
dim(vacc_noncov) <- c(vnc_row, 2)
dim(vacc_child_age) <- vnc_row
dim(vcu_noncov) <- c(vnc_row, 2)
dim(vacc_cu_age) <- vnc_row

## extra outputs
output(time) <- TRUE
output(Ntotal[,,]) <- TRUE
output(inf_1[,,]) <- TRUE
output(FOI1) <- TRUE
output(R0t_1) <- TRUE
output(Delta[]) <- TRUE
output(Kc[]) <- TRUE
output(eip[]) <- TRUE
