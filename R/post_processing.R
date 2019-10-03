
#------------------------------------------------

#' post_processing

#'
#' \code{post_processing} post processes model outputs and saves plots of diagnostics
#'
#' @param dat list of model outputs from the model run.
#'
#' @inheritParams save_plot
#'
#' @importFrom stats setNames
#'
#' @importFrom rlang .data
#'
#' @export


post_processing <- function(dat) {

  weekly_lag_time <- 14

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

  new_levs <- c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Sp", "I1p", "R1p", "Weekly infections/1000")

  levels(df_H_melt$diagnostic) <- new_levs

  df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = c("Susceptibles", "Infectious", "Recovered", "Total population", "Births", "Incidence of infections", "Cumulative incidence", "Weekly infections/1000", "Sp", "I1p", "R1p"))

  diagno_1 <- subset(df_H_melt, diagnostic %in% c("Sp", "I1p", "R1p"))

  diagno_2 <- subset(df_H_melt, diagnostic %in% c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Weekly infections/1000"))

  list("compartments" = droplevels(diagno_1), "demographics" = droplevels(diagno_2))

}
