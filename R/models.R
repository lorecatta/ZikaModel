
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib ZikaModel, .registration = TRUE
#' @import ggplot2
#' @import coda
#' @import scales
#' @import patchwork
#' @import utils
#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("ZikaModel", libpath)
}

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
