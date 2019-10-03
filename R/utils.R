# -----------------------------------------------------------------------------

# Calculates lagged differences using \code{diff}

#' \code{lag_diff} returns a vector of lagged difference values with initial padding.
#'
#' @param x vector of numeric values
#'
#' @param lag the lag value (numeric)
#'
#' @export


lag_diff <- function(x, lag) {

  ret1 <- diff(x, lag = lag, differences = 1)

  pad <- rep(0, lag)

  c(pad, ret1)

}
