
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# -----------------------------------------------------------------------------

## Index locations of outputs in odin model

#' @noRd
odin_index <- function(model) {
  n_out <- environment(model$initialize)$private$n_out %||% 0
  n_state <- length(model$initial(0))
  model$transform_variables(seq_len(1L + n_state + n_out))
}

# -----------------------------------------------------------------------------

## get one odin output

#' @noRd
unpack_odin <- function(x, var_to_unpack) {
  index <- odin_index(x$model)
  temp <- x$output[,unlist(index[var_to_unpack])]
  array(temp, dim = c(dim(temp)[1], x$parameters$na, 2, x$parameters$NP))
}
