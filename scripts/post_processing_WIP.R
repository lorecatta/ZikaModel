
# define parameters -----------------------------------------------------------


out_dir <- file.path("figures", "deterministic")

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

diagno_hum_summary <- c("R0" = "R0t_1av",
                        "FOI on humans" = "FOI1av",
                        "Averaged mortality rate per patch" = "Deltaav",
                        "Mean daily mortality rate" = "DeltaMean")

extra_diagno_hum <- c("Weekly infections/1000" = "w_Ir_infections")

diagno_mos <- c("Susceptibles wt" = "Mwt_S",
                "Infected 1 wt" = "Mwt_E1",
                "Infected 2 wt" = "Mwt_E2",
                "Infectious wt" = "Mwt_I1",
                "Total wt" = "Mwt_tot",
                "Larvae wt birth" = "Lwt_birth",
                "Larvae wt" = "Lwt",
                "Susceptibles wb" = "Mwb_S",
                "Infected 1 wb" = "Mwb_E1",
                "Infected 2 wb" = "Mwb_E2",
                "Infectious wb" = "Mwb_I1",
                "Total wb" = "Mwb_tot",
                "Larvae wb birth" = "Lwb_birth",
                "Larvae wb" = "Lwb")

diagno_mos_summary <- c("FOI on wt mos" = "Mwt_FOI1av")

mos_comparts <- c("Mwt_S",
                  "Mwt_E1",
                  "Mwt_E2",
                  "Mwt_I1",
                  "Mwt_tot",
                  "Mwb_S",
                  "Mwb_E1",
                  "Mwb_E2",
                  "Mwb_I1",
                  "Mwb_tot")


# post processing -------------------------------------------------------------


out <- model_run$dat

tt <- seq(1, dim(out$S)[1], 1)

mos <- setNames(out[diagno_mos[diagno_mos %in% mos_comparts]], diagno_mos[diagno_mos %in% mos_comparts])

dia_hum <- setNames(out[c(diagno_hum, diagno_hum_summary)],
                    c(diagno_hum, diagno_hum_summary))

dia_mos <- setNames(out[c(diagno_mos, diagno_mos_summary)],
                    c(diagno_mos, diagno_mos_summary))

# sum across patches
mossum <- lapply(mos, function(x){apply(x, sum, MARGIN = 1)})

dia_humsum <- lapply(dia_hum[diagno_hum],
                     function(x){apply(x, sum, MARGIN = 1)})

dia_humsum <- c(dia_humsum, dia_hum[diagno_hum_summary])

dia_mossum <- lapply(dia_mos[diagno_mos],
                     function(x) {apply(x, sum, MARGIN = 1)})

dia_mossum <- c(dia_mossum, dia_mos[diagno_mos_summary])


# plot wild type mosquitoes compartments --------------------------------------


mat_Mwt <- do.call("cbind", mossum[mos_comparts[1:4]])

mat_Mwt <- mat_Mwt / mossum$Mwt_tot

mat_Mwt[is.na(mat_Mwt)] <- 0

df_Mwt <- as.data.frame(mat_Mwt)

df_Mwt$time <- tt

df_Mwt_melt <- melt(df_Mwt,
                    id.vars = "time",
                    variable.name = "compartment")

wt_mos_comp_plot <- plot_compartments(df_Mwt_melt,
                                      names(diagno_mos[1:4]),
                                      "SEIR Zika model - wild type mosquitoes")

save_plot(wt_mos_comp_plot,
          "figures",
          "compartments_mosquitoes_wt",
          wdt = 17,
          hgt = 12)


# plot wolbachia mosquitoes compartments --------------------------------------


mat_Mwb <- do.call("cbind", mossum[mos_comparts[6:9]])

mat_Mwb <- mat_Mwb / mossum$Mwb_tot

mat_Mwb[is.na(mat_Mwb)] <- 0

df_Mwb <- as.data.frame(mat_Mwb)

df_Mwb$time <- tt

df_Mwb_melt <- melt(df_Mwb,
                    id.vars = "time",
                    variable.name = "compartment")

wb_mos_comp_plot <- plot_compartments(df_Mwb_melt,
                                      names(diagno_mos[8:11]),
                                      "SEIR Zika model - wolbachia mosquitoes")

save_plot(wb_mos_comp_plot,
          "figures",
          "compartments_mosquitoes_wb",
          wdt = 17,
          hgt = 12)


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
                 out_dir,
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
                 out_dir,
                 "diagnostics_mosquitoes",
                 names(c(diagno_mos, diagno_mos_summary)),
                 no_pages = 2)


# plot by patch ---------------------------------------------------------------

## Proportion of Susceptibles

prop_Sp_df <- as.data.frame(out$prop_Sp)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "Proportion of Susceptibles per patch") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "propSp",
          wdt = 17,
          hgt = 9)

## FOI1p

prop_Sp_df <- as.data.frame(out$FOI1p)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "FOI1p") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "FOI1p",
          wdt = 17,
          hgt = 9)

## Proportion of Infected

prop_Sp_df <- as.data.frame(out$Y1T)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "Proportion of Infected per patch") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "propInfected",
          wdt = 17,
          hgt = 9)

## R0t_1

prop_Sp_df <- as.data.frame(out$R0t_1)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "R0t_1") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "R0t_1",
          wdt = 17,
          hgt = 9)

## Mwt_tot

prop_Sp_df <- as.data.frame(out$Mwt_tot)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "Mwt_tot") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "Mwt_tot",
          wdt = 17,
          hgt = 9)


## M_tot

prop_Sp_df <- as.data.frame(out$M_tot)
colnames(prop_Sp_df) <- seq_len(21)
prop_Sp_df$time <- tt
prop_Sp_df_melt <- melt(prop_Sp_df,
                        id.vars = "time",
                        variable.name = "Patch")
brks <- seq(from = 1, to = max(tt), by = 364 * 10)
prop_Sp_plot <- ggplot(prop_Sp_df_melt, aes(x = time, y = value, colour = Patch)) +
  geom_line(size = 0.4) +
  scale_y_continuous(name = "M_tot") +
  scale_x_continuous(name = "Years", breaks = brks, labels = round(brks / 364))

save_plot(prop_Sp_plot,
          out_dir,
          "M_tot",
          wdt = 17,
          hgt = 9)

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
#                         out_dir,
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
# plot_compartments(df_melt_2, out_dir, out_fig_name)
#
