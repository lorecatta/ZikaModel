
library(stochZika)


# define parameters -----------------------------------------------------------


out_fig_path <- file.path("figures", "zika_model")

diagno_hum <- c("Susceptibles" = "S",
                "Infectious" = "I1",
                "Recovered" = "R1",
                "Total population" = "Ntotal",
                "Births" = "births",
                "Deaths" = "deaths",
                "Incidence of infections" = "inf_1",
                "Cumulative infections" = "sinf1",
                "Fraction of infected h per patch" = "Y1T",
                "Fraction of infectious h per patch" = "infectious1",
                "Incidence of cases" = "disease1",
                "Cumulative cases" = "disease1inc")
                # "Suscep Out" = "O_S",
                # "Suscep ageing" = "age_S",
                # "Infectious Out" = "O_I1",
                # "Infectious recovered" = "recov1",
                # "Infectious ageing" = "age_I1",
                # "Recovered Out" = "O_R1",
                # "Recovered ageing" = "age_R1")

diagno_hum_summary <- c("R0" = "R0t_1av",
                        "FOI on humans" = "FOI1av")
                        #"Averaged mortality rate per patch" = "Deltaav",
                        #"Mean daily mortality rate" = "DeltaMean")

extra_diagno_hum <- c("Weekly infections/1000" = "w_Ir_infections")

diagno_mos <- c("Susceptibles wt" = "Mwt_S",
                "Infected 1 wt" = "Mwt_E1",
                "Infected 2 wt" = "Mwt_E2",
                "Infectious wt" = "Mwt_I1",
                "Total wt" = "Mwt_tot",
                "Larvae wt birth" = "Lwt_birth",
                "Larvae wt" = "Lwt",
                # "Susceptibles wb" = "Mwb_S",
                # "Infected 1 wb" = "Mwb_E1",
                # "Infected 2 wb" = "Mwb_E2",
                # "Infectious wb" = "Mwb_I1",
                # "Total wb" = "Mwb_tot",
                # "Larvae wb birth" = "Lwb_birth",
                # "Larvae wb" = "Lwb",
                "Incidence of infected wt mos" = "Mwt_inf1",
                # "Incidence of infected 1 wt mos leaving" = "O_Mwt_E1",
                "Incidence of suscept wt mos leaving" = "O_Mwt_S",
                "Incidence of wt larvae maturing to adults" = "Lwt_mature")

diagno_mos_summary <- c("FOI on wt mos" = "Mwt_FOI1av")

mos_comparts <- c("Mwt_S",
                  "Mwt_E1",
                  "Mwt_E2",
                  "Mwt_I1",
                  "Mwt_tot")
                  # "Mwb_S",
                  # "Mwb_E1",
                  # "Mwb_E2",
                  # "Mwb_I1",
                  # "Mwb_tot")


# odin ------------------------------------------------------------------------


out <- model_run$


# post processing -------------------------------------------------------------

# diagnostics
hum <- setNames(out[diagno_hum[1:5]], diagno_hum[1:5])

mos <- setNames(out[diagno_mos[diagno_mos %in% mos_comparts]], diagno_mos[diagno_mos %in% mos_comparts])

dia_hum <- setNames(out[c(diagno_hum, diagno_hum_summary)],
                    c(diagno_hum, diagno_hum_summary))

dia_mos <- setNames(out[c(diagno_mos, diagno_mos_summary)],
                    c(diagno_mos, diagno_mos_summary))

# sum across age groups, patches and vaccine statuses
humsum <- lapply(hum, function(x){apply(x, sum, MARGIN = 1)})

# sum across patches
mossum <- lapply(mos, function(x){apply(x, sum, MARGIN = 1)})

dia_humsum <- lapply(dia_hum[diagno_hum],
                     function(x){apply(x, sum, MARGIN = 1)})

dia_humsum <- c(dia_humsum, dia_hum[diagno_hum_summary])

dia_mossum <- lapply(dia_mos[diagno_mos],
                     function(x) {apply(x, sum, MARGIN = 1)})

dia_mossum <- c(dia_mossum, dia_mos[diagno_mos_summary])


# plot human compartments -----------------------------------------------------


mat_H <- do.call("cbind", humsum[c("S", "I1", "R1")])

mat_H <- mat_H / humsum$Ntotal

mat_H[is.na(mat_H)] <- 0

df_H <- as.data.frame(mat_H)

df_H$time <- tt

df_H_melt <- melt(df_H,
                  id.vars = "time",
                  variable.name = "compartment")

plot_compartments(df_H_melt,
                  out_fig_path,
                  "compartments_human.png",
                  names(diagno_hum[1:4]),
                  "SEIR Zika model - human states")


# plot wild type mosquitoes compartments --------------------------------------


mat_Mwt <- do.call("cbind", mossum[mos_comparts[1:4]])

mat_Mwt <- mat_Mwt / mossum$Mwt_tot

mat_Mwt[is.na(mat_Mwt)] <- 0

df_Mwt <- as.data.frame(mat_Mwt)

df_Mwt$time <- tt

df_Mwt_melt <- melt(df_Mwt,
                    id.vars = "time",
                    variable.name = "compartment")

plot_compartments(df_Mwt_melt,
                  out_fig_path,
                  "compartments_mosquitoes_wt.png",
                  names(diagno_mos[1:4]),
                  "SEIR Zika model - wild type mosquitoes")


# plot wolbachia mosquitoes compartments --------------------------------------


# mat_Mwb <- do.call("cbind", mossum[mos_comparts[6:9]])
#
# mat_Mwb <- mat_Mwb / mossum$Mwb_tot
#
# mat_Mwb[is.na(mat_Mwb)] <- 0
#
# df_Mwb <- as.data.frame(mat_Mwb)
#
# df_Mwb$time <- tt
#
# df_Mwb_melt <- melt(df_Mwb,
#                     id.vars = "time",
#                     variable.name = "compartment")
#
# plot_compartments(df_Mwb_melt,
#                   out_fig_path,
#                   "compartments_mosquitoes_wb.png",
#                   names(diagno_mos[5:8]),
#                   "SEIR Zika model - wolbachia mosquitoes")


# plot human diagnostics ------------------------------------------------------


mat_diagnostics <- do.call("cbind", dia_humsum)

df_diagnostics <- as.data.frame(mat_diagnostics)

df_diagnostics$time <- tt

df_diagnostics$w_Ir_infections <- c(rep(0, 7),
                                    diff(df_diagnostics$sinf1,
                                         lag = 7,
                                         differences = 1))

df_diagnostics$w_Ir_infections <- df_diagnostics$w_Ir_infections / df_diagnostics$Ntotal * 1000

# new_array <- array(0, dim = dim(out$sinf1))
#
# for (i in seq_len(11)){
#
#   for (j in seq_len(2)){
#
#     for (k in seq_len(21)){
#
#       new_array[,i,j,k] <- c(rep(0, 7),
#                              diff(out$sinf1[,i,j,k], lag = 7, differences = 1))
#
#     }
#   }
# }
#
# new_array_2 <- array(0, dim = dim(out$sinf1))
#
# new_array_2[out$Ntotal != 0] <- (new_array[out$Ntotal != 0] / out$Ntotal[out$Ntotal != 0]) * 1000
#
# df_diagnostics$w_Ir_infections <- apply(new_array_2, sum, MARGIN = 1)

df_diagnostics_melt <- melt(df_diagnostics,
                            id.vars = "time",
                            variable.name = "diagnostics")

plot_diagnostics(df_diagnostics_melt,
                 out_fig_path,
                 "diagnostics_humans",
                 names(c(diagno_hum, diagno_hum_summary, extra_diagno_hum)),
                 no_pages = 2)


# plot mosquito diagnostics ---------------------------------------------------


mat_diagnostics_mos <- do.call("cbind", dia_mossum)

df_diagnostics_mos <- as.data.frame(mat_diagnostics_mos)

df_diagnostics_mos$time <- tt

df_diagnostics_mos_melt <- melt(df_diagnostics_mos,
                                id.vars = "time",
                                variable.name = "diagnostics")

plot_diagnostics(df_diagnostics_mos_melt,
                 out_fig_path,
                 "diagnostics_mosquitoes",
                 names(c(diagno_mos, diagno_mos_summary)),
                 no_pages = 2)


# plot diagnostics by age groups ----------------------------------------------


# hum_by_age <- setNames(out[c("deathrt", "agert")], c("mortality rate", "ageing rate"))
#
# df_hum_by_age <- lapply(hum_by_age, function(x) as.data.frame(x))
#
# df_hum_by_age <- lapply(df_hum_by_age, give_col_names)
#
# df_hum_by_age <- lapply(df_hum_by_age, add_time_var)
#
# df_diagnostics_melt_1 <- lapply(df_hum_by_age,
#                                 melt,
#                                 id.vars = "time",
#                                 variable.name = "age_group")
#
# df_diagnostics_melt_2 <- lapply(seq_along(df_diagnostics_melt_1),
#                                 add_diagno_name_var,
#                                 b = df_diagnostics_melt_1)
#
# df_diagnostics_by_age <- do.call("rbind", df_diagnostics_melt_2)
#
# df_diagnostics_by_age$diagnostics <- factor(df_diagnostics_by_age$diagnostics,
#                                             levels = c("mortality rate", "ageing rate"))
#
# plot_diagnostics_by_age(df_diagnostics_by_age,
#                         out_fig_path,
#                         "diagnostics_human_by_age",
#                         c("mortality rate", "ageing rate"))
#
# ggplot(df_diagnostics_by_age, aes(x = time, y = value, colour = age_group)) +
#   geom_line(size = 0.4) +
#   facet_wrap(~ diagnostics,
#              ncol = 1,
#              nrow = 2,
#              scales = "free_y")


# # sum across the fourth dimension - patches
# p_humsum <- lapply(hum, function(x){apply(x, sum, MARGIN = 1:3)})
#
# # sum across the second dimension - age groups
# a_p_humsum <- lapply(p_humsum, function(x){apply(x, sum, MARGIN = c(1,3))})
#
# mat <- do.call("cbind", a_p_humsum[1:4])
#
# mat[, c(1,3,5,7)] <- mat[, c(1,3,5,7)] / a_p_humsum$Ntot[, 1]
#
# mat[, c(2,4,6,8)] <- mat[, c(2,4,6,8)] / a_p_humsum$Ntot[, 2]
#
# mat[is.na(mat)] <- 0
#
# df <- as.data.frame(mat)
#
# all_combs <- expand.grid(c("nv", "v"), comp_names)
#
# col_names <- paste(all_combs$Var2, all_combs$Var1, sep = "_")
#
# colnames(df) <- col_names
#
# df$time <- tt
#
# df_melt <- reshape(data = df,
#                    idvar = "time",
#                    varying = col_names,
#                    timevar = "vaccine",
#                    sep = "_",
#                    direction = "long")
#
# df_melt_2 <- melt(df_melt,
#                   id.vars = c("time", "vaccine"),
#                   variable.name = "compartment")
#
# plot_compartments(df_melt_2, out_fig_path, out_fig_name)
#
