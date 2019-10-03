function () {

  tt <- out$TIME
  time <- max(tt)

  diagno_hum <- c("S", "I1", "R1", "births", "inf_1")

  dia_hum <- setNames(out[diagno_hum], diagno_hum)

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

  df_H <- as.data.frame(mat_H)

  # rate of total weekly infections
  df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, 14)

  df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

  df_H$time <- tt
  df_H_melt <- melt(df_H,
                    id.vars = "time",
                    variable.name = "diagnostic")

  diagno_levs <- c("S", "I1", "R1", "Nt", "births", "inf_1", "inf_1_cum", "wIR_inf")

  df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = diagno_levs, labels = diagno_levs)

  ret <- plot_diagnostics(df_H_melt,
                          "human_diagnostics",
                          diagno_levs)

  df_prop <- as.data.frame(prop)

  df_prop$time <- tt
  df_prop_melt <- melt(df_prop,
                       id.vars = "time",
                       variable.name = "compartment")

  diagno_levs_2 <- c("Sp", "I1p", "R1p")

  df_prop_melt$compartment <- factor(df_prop_melt$compartment, levels = diagno_levs_2, labels = diagno_levs_2)

  ret <- plot_diagnostics(df_H_melt,
                          "human_diagnostics",
                          diagno_levs)

  ret2 <- plot_compartments(df_prop_melt,
                            c("Susceptibles", "Infectious", "Recovered"))

  list("diagnostics" = ret, "proportions" = ret2, "dat" = out)

}
