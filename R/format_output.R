
# -----------------------------------------------------------------------------

#' The function formats human-realetd model outputs as a long data frame.
#'
#' @title Format human model outputs
#'
#' @param x Zika_model_simulation object.
#'
#' @param var_select Vector of human compartment names, \code{c("S", "I1", "R1")}. In
#'   addition a number of additional variables can be requested. These include:
#' \itemize{
#'       \item{"Ntotal"}{ Total number of people }
#'       \item{"inf_1"}{ Daily infections }
#'       \item{"MC"}{ Daily microcephaly cases }
#'       \item{"inf_1_w"}{ Weekly infections per 1000 }
#'       \item{"MC_w"}{ Weekly microcephaly cases per 1000 }
#'       }
#'
#' @param keep name of variable to stratify by
#'   (allowed are \code{c("patch", "vaccine", "all")}. Default is no stratification)
#'
#' @param patch_id ID of patch to keep in the summary. Default to keep all patches.
#'
#' @return Formatted long data.frame
#'
#' @export
format_output_H <- function(x,
                            var_select = NULL,
                            keep = NULL,
                            patch_id = NULL) {

  if(!is.null(keep) && keep == "all") {

    if((!is.null(var_select) && (length(var_select) > 1)) | is.null(var_select))
      stop("Can not select more than one compartment when stratifying by both patch and vaccine status")
  }

  # Get model details
  nt <- nrow(x$output)
  index <- odin_index(x$model)

  # Extracting relevant columns for compartment variables
  # -> if var_select = NULL extract all compartments
  # -> if var_select = names extract specific compartments or summary variables

  # Summarise
  # -> if keep = NULL calculate totals across age, vaccine status and patches
  # -> if keep = patch summarise by patch
  # -> if keep = vaccine summarise by vaccine status

  compartments <- c("S", "I1", "R1")

  summary_vars_present <- c("Ntotal", "inf_1")

  summary_vars_not_present <- c("MC")

  weekly_summary_vars <- c("inf_1_w", "MC_w")

  compartments_requested <- var_select[var_select %in% compartments]
  compartments_requested <- if (identical(compartments_requested, character(0))) NULL else compartments_requested

  summary_vars_present_requested <- var_select[var_select %in% summary_vars_present]
  summary_vars_present_requested <- if (identical(summary_vars_present_requested, character(0))) NULL else summary_vars_present_requested

  summary_vars_not_present_requested <- var_select[var_select %in% summary_vars_not_present]
  summary_vars_not_present_requested <- if (identical(summary_vars_not_present_requested, character(0))) NULL else summary_vars_not_present_requested

  weekly_summary_vars_requested <- var_select[var_select %in% weekly_summary_vars]
  weekly_summary_vars_requested <- if (identical(weekly_summary_vars_requested, character(0))) NULL else weekly_summary_vars_requested

  requested_vars <- c(compartments_requested,
                      summary_vars_present_requested,
                      summary_vars_not_present_requested,
                      weekly_summary_vars_requested)

  if (is.null(var_select)) {

    # here the option of summarising by patch or vaccine status
    compartments_output_list <- lapply(compartments, function(j) {

      temp <- x$output[,unlist(index[j])]
      temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$na, 2, x$parameters$NP))
      sum_across_array_dims(temp_array, keep)

    })

    names(compartments_output_list) <- compartments

    output_list <- compartments_output_list

  } else if (!is.null(requested_vars)) {

    output_list <- list()

    vars_present <- c(compartments_requested, summary_vars_present_requested)

    if (!is.null(vars_present)) {

      compartments_output_list <- lapply(vars_present, function(j) {

        temp <- x$output[,unlist(index[j])]
        temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$na, 2, x$parameters$NP))
        sum_across_array_dims(temp_array, keep)

      })

      names(compartments_output_list) <- vars_present

      output_list <- c(output_list, compartments_output_list)

    }

    if ("MC" %in% requested_vars) {

      MC <- calculate_microcases(x)

      MC_output <- sum_across_array_dims(MC, keep)

      output_list <- c(output_list, MC = list(MC_output))

    }

    if ("MC_w" %in% requested_vars) {

      MC <- calculate_microcases(x)

      Ntotal <- unpack_odin(x, "Ntotal")

      sum_MC <- sum_across_array_dims(MC, keep)

      sum_Ntotal <- sum_across_array_dims(Ntotal, keep)

      cumsum_sum_MC <- cumsum_across_array_dims(sum_MC, keep)

      MC_w <- calculate_incidence(cumsum_sum_MC, sum_Ntotal, 7)

      output_list <- c(output_list, MC_w = list(MC_w))

    }

    if ("inf_1_w" %in% requested_vars) {

      inf_1 <- unpack_odin(x, "inf_1")

      Ntotal <- unpack_odin(x, "Ntotal")

      sum_inf_1 <- sum_across_array_dims(inf_1, keep)

      sum_Ntotal <- sum_across_array_dims(Ntotal, keep)

      cumsum_sum_inf_1 <- cumsum_across_array_dims(sum_inf_1, keep)

      inf_1_w <- calculate_incidence(cumsum_sum_inf_1, sum_Ntotal, 7)

      output_list <- c(output_list, inf_1_w = list(inf_1_w))

    }

  } else {

    output_list <- list()

  }

  vars <- names(output_list)

  if(is.null(keep)) {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      compartment = as.character(mapply(rep, vars, nt)),
                      y = unlist(output_list, use.names = FALSE))

  } else if (keep == "patch") {

    out <- data.frame(t = rep(as.numeric(x$output[,index$time]), x$parameters$NP),
                      patch = rep(seq_len(x$parameters$NP), each = nt),
                      compartment = as.character(mapply(rep, vars, x$parameters$NP*nt)),
                      y = unlist(output_list, use.names = FALSE))

  } else if (keep == "vaccine") {

    out <- data.frame(t = rep(as.numeric(x$output[,index$time]), 2),
                      vaccine = rep(seq_len(2), each = nt),
                      compartment = as.character(mapply(rep, vars, 2*nt)),
                      y = unlist(output_list, use.names = FALSE))

  } else if (keep == "all") {

    out <- data.frame(t = rep(as.numeric(x$output[,index$time]), x$parameters$na*2*x$parameters$NP),
                      age = as.factor(rep(seq_len(x$parameters$na), each = nt)),
                      vaccine = rep(seq_len(2), each = nt*x$parameters$na),
                      patch = rep(seq_len(x$parameters$NP), each = nt*x$parameters$na*2),
                      compartment = as.character(mapply(rep, vars, nt*x$parameters$na*2*x$parameters$NP)),
                      y = unlist(output_list, use.names = FALSE))

    if(!is.null(patch_id)) {

      out <- subset(out, patch == patch_id)

    }

  }

  out

}


# -----------------------------------------------------------------------------

#' The function formats mosquito-related model outputs as a long data frame.
#'
#' @title Format mosquito model outputs
#'
#' @param x Zika_model_simulation object.
#'
#' @param var_select Vector of mosquito compartment names,
#' e.g. \code{c("Lwt", "Mwt_S", "Mwb_E1")}, or mosquito-related variables,
#' e.g. \code{c("Kc", "Delta", "R0t_1")}.
#'
#' @param keep name of variable to stratify by
#'   (only allowed \code{"patch"}. Default is no stratification)
#'
#' @return Formatted long data.frame
#'
#' @export
format_output_M <- function(x,
                            var_select = NULL,
                            keep = NULL) {

  # Get model details
  nt <- nrow(x$output)
  index <- odin_index(x$model)

  # Extracting relevant columns for compartment variables
  # -> if var_select = NULL extract all compartments
  # -> if var_select = names specific compartments, extract those

  # Summarise
  # -> if keep = NULL calculate totals across patches (sum or mean depending on the variable)
  # -> if keep = patch summarise by patch

  compartments <- c("Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1",
                    "Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1")

  summary_vars <- c("Lwt", "Lwb", "Kc", "eip", "Delta", "R0t_1", "FOI1")

  compartments_requested <- var_select[var_select %in% compartments]
  compartments_requested <- if (identical(compartments_requested, character(0))) NULL else compartments_requested

  summary_vars_requested <- var_select[var_select %in% summary_vars]
  summary_vars_requested <- if (identical(summary_vars_requested, character(0))) NULL else summary_vars_requested

  requested_vars <- c(compartments_requested, summary_vars_requested)

  if (is.null(var_select)) {

    # here the option of summarising by patch
    output_list <- lapply(compartments, function(j) {

      temp <- x$output[,unlist(index[j])]
      temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$NP))
      sum_across_array_dims(temp_array, keep, j)

    })

    names(output_list) <- compartments

  } else if (!is.null(requested_vars)) {

    output_list <- lapply(requested_vars, function(j) {

      temp <- x$output[,unlist(index[j])]
      temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$NP))
      sum_across_array_dims(temp_array, keep, j)

    })

    names(output_list) <- requested_vars

  } else {

    output_list <- list()

  }

  vars <- names(output_list)

  if(is.null(keep)) {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      compartment = as.character(mapply(rep, vars, nt)),
                      y = unlist(output_list, use.names = FALSE))

  } else if (keep == "patch") {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      patch = rep(seq_len(x$parameters$NP), each = nt),
                      compartment = as.character(mapply(rep, vars, x$parameters$NP*nt)),
                      y = unlist(output_list, use.names = FALSE))

  }

  out

}


# -----------------------------------------------------------------------------

#' The function formats mosquito-related proportions as a long data frame.
#'
#' @title Format mosquito proportions
#'
#' @param x Zika_model_simulation object.
#'
#' @param keep name of variable to stratify by
#'   (only allowed \code{"patch"}. Default is no stratification)
#'
#' @return Formatted long data.frame
#'
#' @export
format_output_Mprop <- function(x, keep = NULL) {

  # Get model details
  nt <- nrow(x$output)
  index <- odin_index(x$model)

  # Summarise
  # -> if keep = NULL calculate totals across patches (sum or mean depending on the variable)
  # -> if keep = patch summarise by patch

  Mwt_comp <- c("Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1")
  Mwb_comp <- c("Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1")

  M_comp <- c(Mwt_comp, Mwb_comp)

  # here the option of summarising by patch
  output_list <- lapply(c(Mwt_comp, Mwb_comp), function(j) {

    temp <- x$output[,unlist(index[j])]
    temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$NP))
    sum_across_array_dims(temp_array, keep, j)

  })

  names(output_list) <- M_comp

  M_tot <- Reduce(`+`, output_list)

  Mwb_tot <- Reduce(`+`, output_list[Mwb_comp])

  prop_wb <- Mwb_tot / (M_tot + 1e-10)

  vars <- "prop_wb"

  if(is.null(keep)) {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      compartment = rep(vars, nt),
                      y = unlist(prop_wb, use.names = FALSE),
                      stringsAsFactors = FALSE)

  } else if (keep == "patch") {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      patch = rep(seq_len(x$parameters$NP), each = nt),
                      compartment = rep(vars, x$parameters$NP*nt),
                      y = c(prop_wb),
                      stringsAsFactors = FALSE)

  }

  out

}
