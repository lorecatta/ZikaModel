#------------------------------------------------
#' post_processing
#'
#' \code{post_processing} post processes model outputs and saves plots of diagnostics.
#'
#' @param dat The dataframe with all the data output from the model run.
#'
#' @inheritParams save_plot
#'
#' @export


post_processing <- function(dat, out_pth){

  diagno_hum <- c("Susceptibles" = "S",
                  "Infectious" = "I1",
                  "Recovered" = "R1",
                  "Total population" = "Ntotal",
                  "Births" = "births",
                  "Incidence of infections" = "inf_1",
                  "Cumulative infections" = "sinf1",
                  "Fraction of infected per patch" = "Y1T",
                  "Fraction of infectious h per patch" = "infectious1")

  diagno_hum_summary <- c("Average R0t_1 across patches" = "R0t_1av",
                          "FOI on humans (average across p)" = "FOI1av",
                          "Average mortality rate across patches" = "Deltaav",
                          "Average carrying capacity across patches" = "Kcav",
                          "Average eip across patches" = "eipav")

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

  diagno_mos_summary <- c("FOI on wt mos (average across p)" = "Mwt_FOI1av")

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


  tt <- seq(1, dim(dat$S)[1], 1)

  mos <- setNames(dat[diagno_mos[diagno_mos %in% mos_comparts]], diagno_mos[diagno_mos %in% mos_comparts])

  dia_hum <- setNames(dat[c(diagno_hum, diagno_hum_summary)],
                      c(diagno_hum, diagno_hum_summary))

  dia_mos <- setNames(dat[c(diagno_mos, diagno_mos_summary)],
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
            out_pth,
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
            out_pth,
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

  # new_array <- array(0, dim = dim(dat$sinf1))
  #
  # for (i in seq_len(11)){
  #
  #   for (j in seq_len(2)){
  #
  #     for (k in seq_len(21)){
  #
  #       new_array[,i,j,k] <- c(rep(0, 7),
  #                              diff(dat$sinf1[,i,j,k], lag = 7, differences = 1))
  #
  #     }
  #   }
  # }
  #
  # new_array_2 <- array(0, dim = dim(dat$sinf1))
  #
  # new_array_2[dat$Ntotal != 0] <- (new_array[dat$Ntotal != 0] / dat$Ntotal[dat$Ntotal != 0]) * 1000
  #
  # df_diagnostics$w_Ir_infections <- apply(new_array_2, sum, MARGIN = 1)

  df_diagnostics_melt <- melt(df_diagnostics,
                              id.vars = "time",
                              variable.name = "diagnostics")

  plot_diagnostics(df_diagnostics_melt,
                   out_pth,
                   "diagnostics_humans",
                   names(c(diagno_hum, diagno_hum_summary, extra_diagno_hum)))


  # plot mosquito diagnostics ---------------------------------------------------


  mat_diagnostics_mos <- do.call("cbind", dia_mossum)

  df_diagnostics_mos <- as.data.frame(mat_diagnostics_mos)

  df_diagnostics_mos$time <- tt

  df_diagnostics_mos_melt <- melt(df_diagnostics_mos,
                                  id.vars = "time",
                                  variable.name = "diagnostics")

  plot_diagnostics(df_diagnostics_mos_melt,
                   out_pth,
                   "diagnostics_mosquitoes",
                   names(c(diagno_mos, diagno_mos_summary)))

}
