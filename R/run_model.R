
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
#' @return list of model outputs
#'
#' @export


run_model <- function(agec,
                      death,
                      nn_links,
                      amplitude_phases,
                      time,
                      season = FALSE) {

  mpl <- model_param_list_create(season = season)

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec = agec,
                                        death = death,
                                        nn_links = nn_links,
                                        amplitude_phases = amplitude_phases,
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
  mod$transform_variables(mod_run)

}
