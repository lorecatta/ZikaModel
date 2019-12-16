
# -----------------------------------------------------------------------------

#' The function processes model outputs to calculate metrics of interest.
#'
#' @title Post-process model outputs
#'
#' @param dat list of model outputs from the model run.
#'
#' @inheritParams model_param_list_create
#'
#' @importFrom stats setNames
#'
#' @export


post_processing <- function(dat, DT) {

  weekly_lag_time <- 7 / DT

  tt <- dat$TIME
  time <- max(tt)

  diagno_hum <- c("S", "I1", "R1", "births", "inf_1")

  dia_hum <- setNames(dat[diagno_hum], diagno_hum)

  # cum sum over the first dimension - specify the dims you want to keep
  # no need for aperm reshaping here
  inf_1_cum <- apply(dia_hum$inf_1, c(2, 3, 4), cumsum)

  Nt <- dia_hum$S + dia_hum$I1 + dia_hum$R1

  dia_hum <- c(dia_hum,
               list(Nt = Nt),
               list(inf_1_cum = inf_1_cum))

  humsum <- lapply(dia_hum, function(x){apply(x, 1, sum)})

  mat_H <- do.call("cbind", humsum)

  prop <- mat_H[, c("S", "I1", "R1")] / mat_H[, "Nt"]

  colnames(prop) <- c("Sp", "I1p", "R1p")

  df_H <- as.data.frame(cbind(mat_H, prop))

  # rate of total weekly infections
  df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, weekly_lag_time)

  df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

  df_H$time <- tt
  df_H_melt <- reshape2::melt(df_H,
                              id.vars = "time",
                              variable.name = "diagnostic")

  current_levs <- c("S", "I1", "R1", "births", "inf_1", "Nt", "inf_1_cum", "Sp", "I1p", "R1p", "wIR_inf")

  # rename
  new_levs <- c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Sp", "I1p", "R1p", "Weekly infections/1000")

  levels(df_H_melt$diagnostic) <- new_levs

  # reorder
  df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = c("Susceptibles", "Infectious", "Recovered", "Total population", "Births", "Incidence of infections", "Cumulative incidence", "Weekly infections/1000", "Sp", "I1p", "R1p"))

  diagno_1 <- subset(df_H_melt, df_H_melt$diagnostic %in% c("Sp", "I1p", "R1p"))

  diagno_2 <- subset(df_H_melt, df_H_melt$diagnostic %in% c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Weekly infections/1000"))

  list("compartments" = droplevels(diagno_1), "demographics" = droplevels(diagno_2))

}



# -----------------------------------------------------------------------------

#' The function processes model outputs to calculate metrics of interest for the
#'  vector population.
#'
#' @title Post-process model outputs for mosquitoes
#'
#' @param dat list of model outputs from the model run.
#'
#' @importFrom stats setNames
#'
#' @export


post_processing_mos <- function(dat) {

  diagno_mos_wt <- c("Lwt", "Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1", "Mwt_tot", "Lwt_birth",
                     "Lwt_mature", "Mwt_inf1")

  diagno_mos_wb <- c("Lwb", "Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1", "Mwb_tot", "Lwb_birth",
                     "Lwb_mature", "Mwb_inf1", "Mwb_intro")

  diagno_mos <- c(diagno_mos_wt, diagno_mos_wb)

  dia_mos <- setNames(dat[diagno_mos], diagno_mos)

  mossum <- lapply(dia_mos, function(x){apply(x, 1, sum)})

  mat_M <- do.call("cbind", mossum)

  df_M <- as.data.frame(mat_M)

  tt <- dat$TIME
  time <- max(tt)

  df_M$time <- tt
  df_M_melt <- melt(df_M,
                    id.vars = "time",
                    variable.name = "diagnostic")

  diagno_levs <- c(diagno_mos_wt, diagno_mos_wb)

  df_M_melt$diagnostic <- factor(df_M_melt$diagnostic, levels = diagno_levs, labels = diagno_levs)

  df_M_melt

}
