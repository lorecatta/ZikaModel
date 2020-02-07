
# -----------------------------------------------------------------------------

#' The function returns a list with the equilibrium initialisation state,
#' odin generator function automatically created given the odin model specified,
#' and list of user-defined parameters.
#'
#' @title Create an odin generator function
#'
#' @inheritParams equilibrium_init_create
#'
#' @inheritParams model_param_list_create
#'
#' @param odin_model_path Character path to odin model.
#'
#' @param params List of user-defined parameters.
#'   Use same names as in \code{model_param_list_create()}. Default = NULL.
#'
#' @return list of generator function, initial state, model parameters and generator
#'
#' @export


create_r_model <- function(odin_model_path,
                           agec,
                           death,
                           nn_links,
                           amplitudes_phases,
                           vaccine_age = NULL,
                           params = NULL) {


  # create model param list with default values
  mpl <- create_model_param_list()

  # replace default parameters with user-defined ones
  if(!is.null(params)) {

    mpl[names(params)] <- params

  }

  # generate initial state variables from equilibrium solution
  state_init <- equilibrium_init_create(agec = agec,
                                        death = death,
                                        nn_links = nn_links,
                                        amplitudes_phases = amplitudes_phases,
                                        vaccine_age = vaccine_age,
                                        model_parameter_list = mpl)

  # create odin generator
  gen <- odin::odin(odin_model_path, verbose = FALSE)
  state <- state_init[names(state_init) %in% names(formals(gen))]

  # return mod
  return(list("state" = state, "generator" = gen, "parameters" = mpl))

}
