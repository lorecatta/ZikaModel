# -----------------------------------------------------------------------------

#' The function calculates the mean across patches of a model output,
#' by ignoring the \emph{rest of the world} patch.
#'
#' @title Calculate the mean across patches of a model output
#'
#' @param my_vector vector of model outputs by patch.
#'
#' @export
calculate_mean_of_patch_variables <- function(my_vector) {

  NP <- length(my_vector)

  ret <- (sum(my_vector) - my_vector[NP]) / (NP - 1)

  return(ret)

}


# -----------------------------------------------------------------------------

#' The function applies a sum over the margin of an array.
#'
#' @title Sum over aray margin
#'
#' @param my_array array to sum.
#'
#' @param keep character vector of which variable / dimension to keep.
#'
#' @param compartment character vector of compartment name.
#'
#' @export
sum_across_array_dims <- function(my_array, keep = NULL, compartment = NULL) {

  summary_vars_to_average <- c("Kc", "eip", "Delta", "R0t_1", "FOI1")

  no_array_dims <- length(dim(my_array))

  if (!is.null(keep) && (no_array_dims == 2 &  keep == "vaccine"))
    stop("Can not summarise mosquito-related variables or compartments by vaccine status")

  if (no_array_dims == 2) {

    if(is.null(keep)) {

      if(compartment %in% summary_vars_to_average) {

        ret <- apply(my_array, 1, calculate_mean_of_patch_variables)

      } else {

        ret <- apply(my_array, 1, sum)

      }

    } else if (keep == "patch") {

      ret <- my_array

    }

  } else {

    if (is.null(keep)) {

      ret <- apply(my_array, 1, sum)

    } else if (keep == "patch") {

      ret <- apply(my_array, c(1, 4), sum)

    } else if (keep == "vaccine") {

      ret <- apply(my_array, c(1, 3), sum)

    } else if (keep == "all") {

      ret <- my_array

    }

  }

  ret

}


# -----------------------------------------------------------------------------

#' The function applies a cumulative sum over the margin of an array.
#'
#' @title Cumulative sum over array margin
#'
#' @param my_array array to sum.
#'
#' @param keep character vector of which variable / dimension to keep
#'
#' @export
cumsum_across_array_dims <- function(my_array, keep = NULL) {

  if (is.null(keep)) {

    ret <- cumsum(my_array)

  } else if (keep == "patch" | keep == "vaccine") {

    ret <- apply(my_array, 2, cumsum)

  }

  ret

}
