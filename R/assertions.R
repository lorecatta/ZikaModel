
# -----------------------------------------------------------------------------
# Dimensions


# -----------------------------------------------------------------------------
# x and y are same length
#' @noRd
assert_same_length <- function(x, y,
                               message =  "%s and %s must be the same length",
                               name_x = deparse(substitute(x)),
                               name_y = deparse(substitute(y))) {
  if (length(x) != length(y)) {
    stop(sprintf(message, name_x, name_y), call. = FALSE)
  }
  return(TRUE)
}


# -----------------------------------------------------------------------------


