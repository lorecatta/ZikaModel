
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
#' @export
format_output_H <- function(x,
                            var_select = NULL,
                            keep = NULL) {

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

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      patch = rep(seq_len(x$parameters$NP), each = nt),
                      compartment = as.character(mapply(rep, vars, x$parameters$NP*nt)),
                      y = unlist(output_list, use.names = FALSE))

  } else if (keep == "vaccine") {

    out <- data.frame(t = as.numeric(x$output[,index$time]),
                      vaccine = rep(seq_len(2), each = nt),
                      compartment = as.character(mapply(rep, vars, 2*nt)),
                      y = unlist(output_list, use.names = FALSE))

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

  compartments <- c("Lwt", "Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1",
                    "Lwb", "Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1")

  summary_vars <- c("Kc", "eip", "Delta", "R0t_1", "FOI1")

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

#' The function processes model outputs to calculate metrics of interest.
#'
#' @title Post-process model outputs
#'
#' @param dat list of model outputs from the model run.
#'
#' @inheritParams parameters_deterministic_model
#'
#' @export
post_processing <- function(dat, DT) {

  weekly_lag_time <- 7 / DT

  tt <- dat$TIME
  time <- max(tt)

  diagno_hum <- c("S", "I1", "R1", "births", "inf_1")

  dia_hum <- setNames(dat[diagno_hum], diagno_hum)

  # cum sum over the first dimension - specify the dims you want to keep
  # no need for aperm reshaping here
  inf_1_cum <- apply(dia_hum$inf_1, c(2, 3, 4), cumsum)

  Nt <- dia_hum$S + dia_hum$I1 + dia_hum$R1

  dia_hum <- c(dia_hum,
               list(Nt = Nt),
               list(inf_1_cum = inf_1_cum))

  humsum <- lapply(dia_hum, function(x){apply(x, 1, sum)})

  mat_H <- do.call("cbind", humsum)

  prop <- mat_H[, c("S", "I1", "R1")] / mat_H[, "Nt"]

  colnames(prop) <- c("Sp", "I1p", "R1p")

  df_H <- as.data.frame(cbind(mat_H, prop))

  # rate of total weekly infections
  df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, weekly_lag_time)

  df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

  df_H$time <- tt
  df_H_melt <- reshape2::melt(df_H,
                              id.vars = "time",
                              variable.name = "diagnostic")

  current_levs <- c("S", "I1", "R1", "births", "inf_1", "Nt", "inf_1_cum", "Sp", "I1p", "R1p", "wIR_inf")

  # rename
  new_levs <- c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Sp", "I1p", "R1p", "Weekly infections/1000")

  levels(df_H_melt$diagnostic) <- new_levs

  # reorder
  df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = c("Susceptibles", "Infectious", "Recovered", "Total population", "Births", "Incidence of infections", "Cumulative incidence", "Weekly infections/1000", "Sp", "I1p", "R1p"))

  diagno_1 <- subset(df_H_melt, df_H_melt$diagnostic %in% c("Sp", "I1p", "R1p"))

  diagno_2 <- subset(df_H_melt, df_H_melt$diagnostic %in% c("Susceptibles", "Infectious", "Recovered", "Births", "Incidence of infections", "Total population", "Cumulative incidence", "Weekly infections/1000"))

  list("compartments" = droplevels(diagno_1), "demographics" = droplevels(diagno_2))

}


# -----------------------------------------------------------------------------

#' The function processes model outputs to calculate metrics of interest for the
#'  vector population.
#'
#' @title Post-process model outputs for mosquitoes
#'
#' @param dat list of model outputs from the model run.
#'
#' @export
post_processing_mos <- function(dat) {

  diagno_mos_wt <- c("Lwt", "Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1", "Mwt_tot", "Lwt_birth",
                     "Lwt_mature", "Mwt_inf1")

  diagno_mos_wb <- c("Lwb", "Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1", "Mwb_tot", "Lwb_birth",
                     "Lwb_mature", "Mwb_inf1", "Mwb_intro")

  diagno_mos <- c(diagno_mos_wt, diagno_mos_wb)

  dia_mos <- setNames(dat[diagno_mos], diagno_mos)

  mossum <- lapply(dia_mos, function(x){apply(x, 1, sum)})

  mat_M <- do.call("cbind", mossum)

  df_M <- as.data.frame(mat_M)

  tt <- dat$TIME
  time <- max(tt)

  df_M$time <- tt
  df_M_melt <- melt(df_M,
                    id.vars = "time",
                    variable.name = "diagnostic")

  diagno_levs <- c(diagno_mos_wt, diagno_mos_wb)

  df_M_melt$diagnostic <- factor(df_M_melt$diagnostic, levels = diagno_levs, labels = diagno_levs)

  df_M_melt

}
