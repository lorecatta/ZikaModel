
# -----------------------------------------------------------------------------

# Create generator for model creation

#' \code{create_r_model} returns list with generator function automatically created given the odin model
#' specified.
#'
#' @inheritParams equilibrium_init_create
#'
#' @inheritParams model_param_list_create
#'
#' @param odin_model_path Character path to odin model
#'
#' @return list of generator function, initial state, model parameters and generator
#'
#' @export


create_r_model <- function(odin_model_path,
                           agec,
                           death,
                           nn_links,
                           amplitude_phases,
                           ...){

  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec = agec,
                                        death = death,
                                        nn_links = nn_links,
                                        amplitude_phases = amplitude_phases,
                                        model_parameter_list = mpl)

  # create odin generator
  gen <- odin::odin(odin_model_path, verbose = FALSE)
  state <- state_init[names(state_init) %in% names(formals(gen))]

  # return mod
  return(list("state" = state, "generator" = gen, "parameters" = mpl))

}
