
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

# -----------------------------------------------------------------------------

#' The function returns a vector of lagged difference values using \code{diff}
#' with initial padding
#'
#' @title Calculate lagged differences
#'
#' @param x vector or matrix of numeric values. In the case of a matrix,
#'   the rows indicate the time-varying obeservations.
#'
#' @param lag the lag value (numeric)
#'
#' @export
lag_diff <- function(x, lag) {

  ret1 <- diff(x, lag = lag, differences = 1)

  if(is.matrix(x) | is.vector(x)) {

    if(is.matrix(x)) {

      pad <- matrix(0, nrow = lag, ncol = ncol(x))

      out <- rbind(pad, ret1)

    }

    if(is.vector(x)) {

      pad <- rep(0, lag)

      out <- c(pad, ret1)

    }

  } else {

    stop("x is not a vector or a matrix")

  }

  out

}
