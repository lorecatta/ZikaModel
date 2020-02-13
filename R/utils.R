
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

    if(is.vectro(x)) {

      pad <- rep(0, lag)

      out <- c(pad, ret1)

    }

  } else {

    stop("x is not a vector or a matrix")

  }

  out

}
