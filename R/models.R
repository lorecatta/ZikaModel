
# -----------------------------------------------------------------------------

#' The function returns a list with
#' the odin model,
#' the function for creating the list of model parameters and
#' the function for running the model.
#'
#' @title Deterministic model creation
#'
#' @export
deterministic_model <- function() {

  model_class <- "deterministic_model"

  deterministic_model <- list(odin_model = odin_model_determ,
                              parameter_func = parameters_deterministic_model,
                              run_func = run_deterministic_model)

  class(deterministic_model) <- c("deterministic", "Zika_model")

  explicit_model

}
