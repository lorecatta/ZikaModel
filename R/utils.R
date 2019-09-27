#' @export
#'

lag_diff <- function(x, lag) {

  ret1 <- diff(x, lag = lag, differences = 1)

  pad <- rep(0, lag)

  c(pad, ret1)

}
