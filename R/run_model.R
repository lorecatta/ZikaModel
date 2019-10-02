
#------------------------------------------------------------------------------

# run_model

#' \code{run_model} runs model using declared agec, death and nn_links
#'
#' @inheritParams equilibrium_init_create
#'
#' @inheritParams model_param_list_create
#'
#' @param time time in days.
#'
#' @importFrom reshape2 melt
#'
#' @importFrom odin odin
#'
#' @export


run_model <- function(agec,
                      death,
                      nn_links,
                      time,
                      season = FALSE){

  mpl <- model_param_list_create(season = season)

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec = agec,
                                        death = death,
                                        nn_links = nn_links,
                                        model_parameter_list = mpl)

  # create odin generator
  odin_model_path <- system.file("extdata/odin_model_determ.R", package = "ZikaModel")

  gen <- odin::odin(odin_model_path,verbose = FALSE)

  state_use <- state_init[names(state_init) %in% names(formals(gen))]

  # create model with initial values
  mod <- gen(user = state_use)
  tt <- seq(0, time, 1)

  # run model
  mod_run <- mod$run(tt)

  # shape output
  out <- mod$transform_variables(mod_run)

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

  mat_H <- cbind(mat_H, prop)

  df_H <- as.data.frame(mat_H)

  # rate of total weekly infections
  df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, 14)

  df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

  df_H$time <- tt
  df_H_melt <- melt(df_H,
                    id.vars = "time",
                    variable.name = "diagnostic")

  diagno_levs <- c("S", "I1", "R1", "Nt", "Sp", "I1p", "R1p", "births", "inf_1", "inf_1_cum", "wIR_inf")

  df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = diagno_levs, labels = diagno_levs)

  ret <- plot_diagnostics(df_H_melt,
                          "human_diagnostics",
                          diagno_levs)

  list("plot" = ret, "dat" = out)

}
