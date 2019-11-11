
# -----------------------------------------------------------------------------

#' The function returns a vector of lagged difference values using \code{diff}
#' with initial padding
#'
#' @title Calculate lagged differences
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
