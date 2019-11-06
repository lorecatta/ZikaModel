
DT <- user()
YL <- user()
TIME <- step * DT

NP <- user()

nn[,] <- user()

age_per <- user()
na <- user()

vnc_row <- na + 1 # position 1 is like position 0 in BM

agec[] <- user()
deathrt[] <- user()

TimeMwtControlOn <- user()
TimeMwtControlOff <- user()
propMwtControl <- user()

MwtCont <- if ((TIME >= TimeMwtControlOn * YL) &&
               (TIME < TimeMwtControlOff * YL)) (1 - propMwtControl) else 1

Rm <- user() # mosquito reproduction number
Epsilon <- user() # larvae development rate
Sigma <- user() # larvae mortality rate
Omega <- user() # 1

# rate at which adult female produce larvae. Rm is fixed.
Gamma <- Rm * DeltaMean * (Epsilon + Sigma) / Epsilon

N_eq[] <- user()
Kappa <- user()
inf_per <- user()

Mwt_cv <- user()
Mwt_mean <- user()

ln_sd <- sqrt(log(1 + Mwt_cv * Mwt_cv))
ln_mean <- log(Mwt_mean) - (ln_sd * ln_sd) / 2

dim(ln_pick) <- NP
ln_pick[] <- ln_mean

Mwt[] <- exp(ln_pick[i])


# -----------------------------------------------------------------------------
# Mean value of mosquito mortality, carrying capacity and EIP
# without the effect of seasonality


DeltaBase <- user()

# mean adult mosquito mortality rate, corrected by the effect of control
DeltaMean <- DeltaBase / MwtCont

# mean larval mosquito carrying capacity. DeltaMean is fixed.
Kc_mean <- DeltaMean *
  ((Epsilon * (Gamma - DeltaMean) / (DeltaMean * Sigma) - 1) ^ (-1 / Omega)) / Epsilon

eip_mean <- user()

# calculate R0 at equilibrium ?
Beta_mh_1 <- user()
Beta_hm_1 <- user()


# -----------------------------------------------------------------------------
# Mean value of mosquito mortality, carrying capacity and EIP
# with the effect of seasonality


season_phase[] <- user()
season_amp[] <- user()
pi <- user()
Delta_season <- user()
Kc_season <- user()
eip_season <- user()

Delta[1:(NP-1)] <- DeltaMean /
  (1 + season_amp[i] * Delta_season * cos(2 * pi * (TIME + season_phase[i]) / YL))
Delta[NP] <- DeltaMean

Kc[1:(NP-1)] <- Mwt[i] * Kc_mean *
  (1 + season_amp[i] * Kc_season * cos(2 * pi * (TIME + season_phase[i]) / YL))
Kc[NP] <- MwtCont * Mwt_mean * Kc_mean *
  (1 + season_amp[i] * Kc_season * cos(2 * pi * (TIME + season_phase[i]) / YL))

eip[] <- eip_mean *
  (1 - season_amp[i] * eip_season * cos(2 * pi * (TIME + season_phase[i]) / YL))

Deltaav <- (sum(Delta[]) - Delta[NP]) / (NP - 1)
Kcav <- (sum(Kc[]) - Kc[NP]) / (NP - 1)
eipav <- (sum(eip[]) - eip[NP]) / (NP - 1)


# -----------------------------------------------------------------------------
# Parameters for interventions



Wb_fM <- user()
Wb_relsusc1 <- user()

Wb_introtime[] <- user()
Wb_introduration <- user()
Wb_introlevel <- user()

Delta_wb[] <- Delta[i] / Wb_fM

# vacc_cu_minage <- user()
# vacc_cu_maxage <- user()
# vacc_cu_coverage <- user()
# vacc_cu_time <- user()
# vacc_cu_rndtime <- user()
#
# vacc_child_age <- user()
# vacc_child_coverage <- user()
# vacc_child_starttime <- user()
# vacc_child_stoptime <- user()



# -----------------------------------------------------------------------------
#
# States of wild type mosquitoes
#
# -----------------------------------------------------------------------------



initial(Lwt[]) <- trunc(Mwt[i] * N_eq[i] * Delta[i] / Epsilon)
initial(Mwt_S[]) <- trunc(Mwt[i] * N_eq[i])
initial(Mwt_E1[]) <- 0
initial(Mwt_E2[]) <- 0
initial(Mwt_I1[]) <- 0

Mwt_tot[] <- Mwt_S[i] + Mwt_E1[i] + Mwt_E2[i] + Mwt_I1[i]

Wb_cyto <- user()
Wb_fF <- user()
Wb_mat <- user()
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

Mwt_FOI1av <- (sum(Mwt_FOI1[]) - Mwt_FOI1[NP]) / (NP - 1)

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



initial(Lwb[]) <- 0
initial(Mwb_S[]) <- 0
initial(Mwb_E1[]) <- 0
initial(Mwb_E2[]) <- 0
initial(Mwb_I1[]) <- 0

Mwb_tot[] <- Mwb_S[i] + Mwb_E1[i] + Mwb_E2[i] + Mwb_I1[i]

dim(M_tot) <- NP
dim(prop_wb) <- NP
M_tot[] <- Mwt_tot[i] + Mwb_tot[i]
prop_wb[] <- Mwb_tot[i] / (M_tot[i] + 1e-10)

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

Mwb_FOI1av <- (sum(Mwb_FOI1[]) - Mwb_FOI1[NP]) / (NP - 1)

O_Mwb_S_prob[] <- max(min(DT * Delta_wb[i] + Mwb_FOI1[i], 1), 0)
O_Mwb_S[] <- Mwb_S[i] * O_Mwb_S_prob[i]

Mwb_inf1_prob[] <- max(min(Mwb_FOI1[i] / (DT * Delta_wb[i] + Mwb_FOI1[i]), 1), 0)
Mwb_inf1[] <- O_Mwb_S[i] * Mwb_inf1_prob[i]

Wb_introrate[] <- Wb_introlevel * Mwt[i] * N_eq[i] * DT / Wb_introduration

Mwb_intro[] <- if((TIME >= Wb_introtime[i] * YL) &&
                  (TIME < Wb_introtime[i] * YL + Wb_introduration)) Wb_introrate[i] else 0

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

Mwt_propinf[] <- (Mwt_E1[i] + Mwt_E2[i] + Mwt_I1[i]) / (Mwt_tot[i] + 1e-10)
Mwb_propinf[] <- (Mwb_E1[i] + Mwb_E2[i] + Mwb_I1[i]) / (Mwb_tot[i] + 1e-10)
M_propinf[] <- ((Mwt_tot[i] + 1e-10) * Mwt_propinf[i] + (Mwb_tot[i] + 1e-10) * Mwb_propinf[i]) / (M_tot[i]+1e-10)



# -----------------------------------------------------------------------------
#
# Human states
#
# -----------------------------------------------------------------------------



init_S[,,] <- user()
initial(S[,,]) <- init_S[i,j,k]

init_I1[,,] <- user()
initial(I1[,,]) <- init_I1[i,j,k]

init_R1[,,] <- user()
initial(R1[,,]) <- init_R1[i,j,k]

Ntotal[,,] <- S[i,j,k] + I1[i,j,k] + R1[i,j,k]

Ntotal_sum[1,1:NP] <- Ntotal[i,1,j] + Ntotal[i,2,j]
Ntotal_sum[2:na,1:NP] <- Ntotal_sum[i-1,j] + Ntotal[i,1,j] + Ntotal[i,2,j]
NTp[] <- Ntotal_sum[na,i]

sum_S[1,1:NP] <- S[i,1,j]
sum_S[2:na,1:NP] <- S[i,1,j] + sum_S[i-1,j]
prop_Sp[] <- sum_S[na,i] / NTp[i]

phi1[] <- user()
Y1[,,] <- phi1[j] * inf_1[i,j,k]

Y1T_sum[1,1:NP] <- Y1[i,1,j] + Y1[i,2,j]
Y1T_sum[2:na,1:NP] <- Y1T_sum[i-1,j] + Y1[i,1,j] + Y1[i,2,j]

Y1T[] <- Y1T_sum[na,i] / NTp[i]



# -----------------------------------------------------------------------------
#
# Keep track of infectious humans
#
# -----------------------------------------------------------------------------



# 2 human incubation and infectious classes (instead of 1)
# are needed for a gamma rather than
# an exponentially distributed human infectious period

initial(incubA[]) <- 0
initial(incubB[]) <- 0
initial(infectiousA[]) <- 0
initial(infectiousB[]) <- 0

update(incubA[]) <- incubA[i] + Y1T[i] - DT * 2 * incubA[i] / incub

update(incubB[]) <- incubB[i] + DT * 2 * (incubA[i] - incubB[i]) / incub

update(infectiousA[]) <- infectiousA[i] +
  DT * 2 * (incubB[i] / incub - infectiousA[i] / inf_per)

update(infectiousB[]) <- infectiousB[i] +
  DT * 2 * (infectiousA[i] - infectiousB[i]) / inf_per

infectious1[] <- infectiousA[i] + infectiousB[i]

incub <- user()

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



Wb_relinf1 <- user()

FOI1p[] <- DT * Beta_mh_1 * Kappa * (Mwt_I1[i] + Mwb_I1[i] * Wb_relinf1) / NTp[i]
FOI1Y[] <- YL * FOI1p[i] / DT

FOI1av <- (sum(FOI1p[]) - FOI1p[NP]) / (NP - 1)

FOI1nn[1:(NP-1)] <- (FOI1p[nn[i,1]] + FOI1p[nn[i,2]] + FOI1p[nn[i,3]] +
                     FOI1p[nn[i,4]] + FOI1p[nn[i,5]] + FOI1p[nn[i,6]] +
                     FOI1p[nn[i,7]] + FOI1p[nn[i,8]]) / 8

propTransGlobal <- user()
propTransNN <- user()
BG_FOI <- user()
FOI1[1:(NP-1)] <- propTransNN * FOI1nn[i] + propTransGlobal * FOI1av +
  (1 - propTransGlobal - propTransNN) * FOI1p[i] + DT * BG_FOI / YL
# FOI1[NP] <- pTG_bigpatch * FOI1av + (1 - pTG_bigpatch) * FOI1p[i]
FOI1[NP] <- 0



# -----------------------------------------------------------------------------
#
# Human compartmental dynamics
#
# -----------------------------------------------------------------------------



agerts <- if (trunc(TIME / age_per) == TIME / age_per) age_per / YL else 0
# agerts=DT/YL
agert[] <- agerts / agec[i] # ageing rate per age group

Nb[] <- user()
births_lambda[] <- DT * Nb[i] / YL
births[] <- births_lambda[i]

rho1[] <- user()
O_S_prob[,,] <- max(min(rho1[j] * FOI1[k] + agert[i] + deathrt[i], 1), 0)
O_S[,,] <- S[i,j,k] * O_S_prob[i,j,k]

inf_1_prob[,,] <- max(min(rho1[j] * FOI1[k] / (rho1[j] * FOI1[k] +
                                                 agert[i] + deathrt[i]), 1), 0)
inf_1[,,] <- O_S[i,j,k] * inf_1_prob[i,j,k]

age_S_trials[,,] <- O_S[i,j,k] - inf_1[i,j,k]
age_S_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_S[,,] <- age_S_trials[i,j,k] * age_S_prob[i]

update(S[1,1,1:NP]) <- trunc(0.5 + births[k] + S[i,j,k] - O_S[i,j,k])
update(S[2:na,1,1:NP]) <- trunc(0.5 + age_S[i-1,j,k] + S[i,j,k] - O_S[i,j,k])

update(S[1,2,1:NP]) <- trunc(0.5 + S[i,j,k] - O_S[i,j,k])
update(S[2:na,2,1:NP]) <- trunc(0.5 + age_S[i-1,j,k] + S[i,j,k] - O_S[i,j,k])

nu <- user()

O_I1_prob[] <- max(min(nu + agert[i] + deathrt[i], 1), 0)
O_I1[,,] <- I1[i,j,k] * O_I1_prob[i]

recov1_prob[] <- max(min(nu / (nu + agert[i] + deathrt[i]), 1), 0)
recov1[,,] <- O_I1[i,j,k] * recov1_prob[i]

age_I1_trials[,,] <- O_I1[i,j,k] - recov1[i,j,k]
age_I1_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_I1[2:vnc_row,1:2,1:NP] <- age_I1_trials[i-1,j,k] * age_I1_prob[i-1]
age_I1[1,1:2,1:NP] <- 0

update(I1[1:na,1:2,1:NP]) <- trunc(0.5 + age_I1[i,j,k] + inf_1[i,j,k] + I1[i,j,k] - O_I1[i,j,k])

O_R1_prob[] <- max(min(agert[i] + deathrt[i], 1), 0)
O_R1[,,] <- R1[i,j,k] * O_R1_prob[i]

age_R1_prob[] <- max(min(agert[i] / (agert[i] + deathrt[i]), 1), 0)
age_R1[2:vnc_row,1:2,1:NP] <- O_R1[i-1,j,k] * age_R1_prob[i-1]
age_R1[1,1:2,1:NP] <- 0

update(R1[1:na,1:2,1:NP]) <- trunc(0.5 + age_R1[i,j,k] + recov1[i,j,k] + R1[i,j,k] - O_R1[i,j,k])



# -----------------------------------------------------------------------------
#
# The outputs you want
#
# -----------------------------------------------------------------------------



# diagnostics for humans
output(agert[]) <- TRUE
output(deathrt[]) <- TRUE
output(births[]) <- TRUE
output(Ntotal[,,]) <- TRUE
output(FOI1av) <- TRUE
output(O_S[,,]) <- TRUE
output(inf_1[,,]) <- TRUE
output(age_S[,,]) <- TRUE
output(O_I1[,,]) <- TRUE
output(recov1[,,]) <- TRUE
output(age_I1[,,]) <- TRUE
output(O_R1[,,]) <- TRUE
output(age_R1[,,]) <- TRUE
output(Y1T[]) <- TRUE
output(Deltaav) <- TRUE
output(Kcav) <- TRUE
output(eipav) <- TRUE
output(FOI1[]) <- TRUE
output(FOI1p[]) <- TRUE
output(FOI1nn[]) <- TRUE
output(prop_Sp[]) <- TRUE
output(FOI1Y[]) <- TRUE
output(inf_1_prob[,,]) <- TRUE
output(Kc[]) <- TRUE
output(eip[]) <- TRUE
output(Delta[]) <- TRUE
output(NTp[]) <- TRUE
output(Mwt_FOI1[]) <- TRUE
output(O_S_prob[,,]) <- TRUE
output(rho1[]) <- TRUE
output(infectious1[]) <- TRUE
output(R0t_1[]) <- TRUE

# diagnostics for mosquitoes
output(Mwt_tot[]) <- TRUE
output(Lwt_birth[]) <- TRUE
output(Lwt_mature[]) <- TRUE
output(Mwt_FOI1av) <- TRUE
output(Mwt_inf1[]) <- TRUE
output(Mwt_propinf) <- TRUE

output(MwtCont) <- TRUE

output(Mwb_tot[]) <- TRUE
output(Lwb_birth[]) <- TRUE
output(Lwb_mature[]) <- TRUE
output(Mwb_FOI1av) <- TRUE
output(Mwb_inf1[]) <- TRUE
output(Mwb_propinf) <- TRUE

output(Wb_introrate[]) <- TRUE
output(Mwb_intro[]) <- TRUE
output(prop_wb[]) <- TRUE

output(M_propinf) <- TRUE

output(TIME) <- TRUE


# -----------------------------------------------------------------------------
#
# define array dimensions (odin requires it)
#
# -----------------------------------------------------------------------------



dim(Mwt) <- NP
dim(season_phase) <- NP
dim(season_amp) <- NP
dim(Delta) <- NP
dim(Kc) <- NP
dim(eip) <- NP
dim(nn) <- c(NP, 8)
dim(N_eq) <- NP
dim(Lwt) <- NP
dim(Mwt_S) <- NP
dim(Mwt_E1) <- NP
dim(Mwt_E2) <- NP
dim(Mwt_I1) <- NP
dim(Mwt_tot) <- NP
dim(Lwt_birth) <- NP
dim(Lwt_birth_lambda) <- NP
dim(L_deathrt) <- NP
dim(L_dr) <- NP
dim(O_Lwt) <- NP
dim(O_Lwt_prob) <- NP
dim(Lwt_mature) <- NP
dim(Lwt_mature_prob) <- NP
dim(Mwt_FOI1) <- NP
dim(O_Mwt_S) <- NP
dim(O_Mwt_S_prob) <- NP
dim(Mwt_inf1) <- NP
dim(Mwt_inf1_prob) <- NP
dim(O_Mwt_E1) <- NP
dim(O_Mwt_E1_prob) <- NP
dim(Mwt_E1_incub) <- NP
dim(Mwt_E1_incub_prob) <- NP
dim(O_Mwt_E2) <- NP
dim(O_Mwt_E2_prob) <- NP
dim(Mwt_E2_incub) <- NP
dim(Mwt_E2_incub_prob) <- NP
dim(O_Mwt_I1) <- NP
dim(O_Mwt_I1_prob) <- NP
dim(init_S) <- c(na, 2, NP)
dim(S) <- c(na, 2, NP)
dim(init_I1) <- c(na, 2, NP)
dim(I1) <- c(na, 2, NP)
dim(init_R1) <- c(na, 2, NP)
dim(R1) <- c(na, 2, NP)
dim(Ntotal) <- c(na, 2, NP)
dim(Ntotal_sum) <- c(na, NP)
dim(NTp) <- NP
dim(sum_S) <- c(na, NP)
dim(prop_Sp) <- NP
dim(Y1) <- c(na, 2, NP)
dim(phi1) <- 2
dim(Y1T_sum) <- c(na,NP)
dim(Y1T) <- NP
dim(infectious1) <- NP
dim(FOI1p) <- NP
dim(FOI1Y) <- NP
dim(FOI1nn) <- NP
dim(FOI1) <- NP
dim(agert) <- na
dim(agec) <- na
dim(Nb) <- NP
dim(births) <- NP
dim(births_lambda) <- NP
dim(O_S) <- c(na, 2, NP)
dim(deathrt) <- na
dim(rho1) <- 2
dim(O_S_prob) <- c(na, 2, NP)
dim(inf_1) <- c(na, 2, NP)
dim(inf_1_prob) <- c(na, 2, NP)
dim(age_S) <- c(na, 2, NP)
dim(age_S_trials) <- c(na, 2, NP)
dim(age_S_prob) <- na
dim(O_I1) <- c(na, 2, NP)
dim(O_I1_prob) <- na
dim(recov1) <- c(na, 2, NP)
dim(recov1_prob) <- na
dim(age_I1) <- c(vnc_row, 2, NP)
dim(age_I1_prob) <- na
dim(age_I1_trials) <- c(na, 2, NP)
dim(O_R1) <- c(na, 2, NP)
dim(O_R1_prob) <- na
dim(age_R1) <- c(vnc_row, 2, NP)
dim(age_R1_prob) <- na
# dim(Y1T_del_inc) <- NP
# dim(Y1T_del_inc_ip) <- NP
dim(incubA) <- NP
dim(incubB) <- NP
dim(infectiousA) <- NP
dim(infectiousB) <- NP

dim(Delta_wb) <- NP
dim(Lwb) <- NP
dim(Mwb_S) <- NP
dim(Mwb_E1) <- NP
dim(Mwb_E2) <- NP
dim(Mwb_I1) <- NP
dim(Mwb_tot) <- NP
dim(R0t_1) <- NP
dim(Lwb_birth) <- NP
dim(Lwb_birth_lambda) <- NP
dim(O_Lwb) <- NP
dim(O_Lwb_prob) <- NP
dim(Lwb_mature) <- NP
dim(Lwb_mature_prob) <- NP
dim(Mwb_FOI1) <- NP
dim(O_Mwb_S) <- NP
dim(O_Mwb_S_prob) <- NP
dim(Mwb_inf1) <- NP
dim(Mwb_inf1_prob) <- NP
dim(Wb_introrate) <- NP
dim(Mwb_intro) <- NP
dim(Wb_introtime) <- NP
dim(O_Mwb_E1) <- NP
dim(O_Mwb_E1_prob) <- NP
dim(Mwb_E1_incub) <- NP
dim(Mwb_E1_incub_prob) <- NP
dim(O_Mwb_E2) <- NP
dim(O_Mwb_E2_prob) <- NP
dim(Mwb_E2_incub) <- NP
dim(Mwb_E2_incub_prob) <- NP
dim(O_Mwb_I1) <- NP
dim(O_Mwb_I1_prob) <- NP
dim(Mwt_propinf) <- NP
dim(Mwb_propinf) <- NP
dim(M_propinf) <- NP
