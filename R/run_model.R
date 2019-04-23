#------------------------------------------------
#' run_model
#'
#' \code{run_model} runs model using declared agec, death and nn_links
#'
#' @inheritParams equilibrium_init_create
#'
#' @param time time in days.
#'   Default=18200
#'
#' @importFrom reshape2 melt
#' @importFrom odin odin
#'
#' @export


run_model <- function(age,
                      death,
                      nn_links,
                      time = 18200){

  mpl <- model_param_list_create()

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec = agec,
                                        death = death,
                                        nn_links = nn_links,
                                        model_parameter_list = mpl)

  # create odin generator
  odin_model_path <- system.file("extdata/odin_model_stoch.R", package = "stochZika")
  gen <- odin::odin(odin_model_path,verbose = FALSE)

  # There are many parameters used that should not be passed through
  # to the model.
  state_use <- state_init[names(state_init) %in% names(formals(gen))]

  # create model with initial values
  mod <- gen(user = state_use)
  tt <- seq(0, time, 1)

  # run model
  mod_run <- mod$run(tt)

  # shape output
  out <- mod$transform_variables(mod_run)
  hum <- list("S" = out$S, "I" = out$I1, "R" = out$R1, "Nt" = out$Ntotal)
  humsum <- lapply(hum, function(x){apply(x, sum, MARGIN = 1)})
  mat_H <- do.call("cbind", humsum[c("S", "I", "R")])
  mat_H <- mat_H / humsum$Nt
  mat_H[is.na(mat_H)] <- 0
  df_H <- as.data.frame(mat_H)
  df_H$time <- tt
  df_H_melt <- melt(df_H,
                    id.vars = "time",
                    variable.name = "compartment")

  ret <- plot_compartments(df_H_melt,
                           c("Susceptibles", "Infectious", "Recovered"),
                           "SEIR Zika model - human states")

  list("plot" = ret, "dat" = out)

}
