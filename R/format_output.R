


# -----------------------------------------------------------------------------

#' The function format model outputs as a long data frame.
#'
#' @title Format model outputs
#'
#' @param x Zika_model_simulation object.
#'
#' @param var_select Vector of compartment names, e.g. \code{c("S", "R")}. In
#'   addition a number of additional variables can be requested. These include:
#' \itemize{
#'       \item{"infections_w"}{ Weekly Infections }
#'       \item{"micro_cases_w"}{ Weekly Microcephaly Cases }
#'       \item{"Kcav"}{ Average mosquito larvae carrying capacity by patch }
#'       \item{"eipav"}{ Average EIP by patch }
#'       \item{"Deltaav"}{ Average adult mosquito mortality rate by patch }
#'       \item{"R0t_1av"}{ Average Rt by patch }
#'       \item{"FOI1av"}{ Average FOI by patch }
#'       }
#' @export
format_output <- function(x,
                          var_select = NULL,
                          keep = NULL) {

  # Get model details
  nt <- nrow(x$output)
  index <- odin_index(x$model)

  browser()

  # Extracting relevant columns for compartment variables
  # -> if var_select = NULL extract all compartments
  # -> if var_select = names specific compartments, extract those

  # Summarise
  # -> if keep = NULL calculate totals across age, vaccine status and patches
  # -> if keep = patch summarise by patch
  # -> if keep = vaccine summarise by vaccine status

  ## calculate number of microcephaly cases
  MC <- calculate_microcases(x)

  patch_vars_to_sum <- c("Lwt", "Mwt_S", "Mwt_E1", "Mwt_E2", "Mwt_I1",
                         "Lwb", "Mwb_S", "Mwb_E1", "Mwb_E2", "Mwb_I1")

  patch_vars_to_average <- c("Kc", "eip", "Delta", "R0t_1", "FOI1")

  human_compartments <- c("S", "I1", "R1")

  inc_vars <- c("infections_w", "micro_cases_w")

  if(is.null(var_select)) {

    patch_vars_to_sum_output_list <- lapply(patch_vars_to_sum, function(j) {

      temp_array <- x$output[,unlist(index[j])]
      apply(temp_array, 1, sum)

    })

    names(patch_vars_to_sum_output_list) <- patch_vars_to_sum

    patch_vars_to_average_output_list <- lapply(patch_vars_to_average, function(j) {

      temp_array <- x$output[,unlist(index[j])]
      apply(temp_array, 1, calculate_mean_of_patch_variables)

    })

    names(patch_vars_to_average_output_list) <- patch_vars_to_average

    # here the option of summarising by patch or vaccine status
    human_compartments_output_list <- lapply(human_compartments, function(j) {

      temp <- x$output[,unlist(index[j])]
      temp_array <- array(temp, dim = c(dim(temp)[1], x$parameters$na, 2, x$parameters$NP))
      sum_across_array_dims(temp_array, keep)

    })

    names(human_compartments_output_list) <- human_compartments

  }

  browser()

  if(is.null(keep)) {

    output_list <- c(human_compartments_output_list)

    mos_out <- data.frame(t = as.numeric(x$output[,index$time]),
                          y = unlist(output_list))

  }


  # Disaggregating var_select into compartments and summary variables
  compartments <- var_select[!(var_select %in% mean_vars)]
  compartments <- if (identical(compartments, character(0))) NULL else compartments
  summaries <- var_select[var_select %in% mean_vars]
  summaries <- if (identical(summaries, character(0))) NULL else summaries


}

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

  # FOI1av <- (sum(FOI1p[]) - FOI1p[NP]) / (NP - 1)
  # R0t_1av <- (sum(R0t_1[]) - R0t_1[NP]) / (NP - 1)
  # Deltaav <- (sum(Delta[]) - Delta[NP]) / (NP - 1)
  # Kcav <- (sum(Kc[]) - Kc[NP]) / (NP - 1)
  # eipav <- (sum(eip[]) - eip[NP]) / (NP - 1)
  # Mwt_FOI1av <- (sum(Mwt_FOI1[]) - Mwt_FOI1[NP]) / (NP - 1)
  # Mwb_FOI1av <- (sum(Mwb_FOI1[]) - Mwb_FOI1[NP]) / (NP - 1)

}

sum_across_array_dims <- function(array_to_sum, keep = NULL) {

  if(is.null(keep)) {

    ret <- apply(array_to_sum, 1, sum)

  } else if(keep == "patch") {

    ret <- apply(array_to_sum, c(1, 4), sum)

  } else if(keep == "vaccine") {

    ret <- apply(array_to_sum, c(1, 3), sum)

  }

  ret

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
