
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

lag_diff_array <- my_fun <- function(my_array, lag_time) {

  dims <- dim(my_array)

  dim_1 <- dims[3]
  dim_2 <- dims[4]

  out <- array(0, dim = dims)

  for (i in seq_len(dim_1)) {

    for (j in seq_len(dim_2)) {

      # browser()

      out[,,i,j] <- lag_diff(my_array[,,i,j], lag_time)

    }

  }

  out

}

calculate_incidence <- function(cum_infections, Ntotal, time_window) {

  n_dims <- length(dim(cum_infections))

  if(time_window == 1) {

    out <- cum_infections

  }

  if(time_window > 1) {

    if(n_dims == 0 | n_dims == 2) {

      out <- lag_diff(cum_infections, time_window)

    }

    if(n_dims == 4) {

      out <- lag_diff_array(cum_infections, time_window)

    }

  }

  ifelse(Ntotal == 0, 0, out / Ntotal * 1000)

}
