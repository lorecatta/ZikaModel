
# -----------------------------------------------------------------------------

#' The function makes a plot of the model compartments.
#'
#' @title Plot human or mosquitoes states
#'
#' @param x Zika_model_simulation object.
#'
#' @param type Character of model compartment type to plot.
#'   allowed are \code{c("H", "M")} for Human and Mosquito.
#'
#' @param var_select Vector of variable names to plot (default is all)
#'
#' @param keep name of variables to stratify by
#'   (allowed are \code{c("patch", "vaccine")} for type H and
#'   \code{"patch"} for type M. Default is no stratification)
#'
#' @export
plot.Zika_model_simulation <- function(x, type, var_select = NULL, keep = NULL) {

  if (type == "H") {

    pds <- format_output_H(x, var_select = var_select, keep = keep)

  }

  if (type == "M") {

    pds <- format_output_M(x, var_select = var_select, keep = keep)

  }

  # Plot
  p <- ggplot2::ggplot()

  if(is.null(keep)) {
    p <- p + ggplot2::geom_line(data = pds,
                                ggplot2::aes(x = .data$t, y = .data$y,
                                             col = .data$compartment))
  } else if (keep == "patch" | keep == "vaccine") {
    p <- p + ggplot2::geom_line(data = pds,
                                ggplot2::aes(x = .data$t, y = .data$y,
                                             col = .data$compartment)) +
      ggplot2::facet_wrap(as.formula(paste("~", keep)),
                          ncol = 4,
                          scales = "free_y")
  }

  # Add remaining formatting
  p <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = FALSE) +
    ggplot2::scale_x_continuous(name = "Time") +
    ggplot2::scale_y_continuous(name = "N", labels = scales::comma) +
    ggplot2::theme_bw()

  return(p)

}


# -----------------------------------------------------------------------------

#' The function makes a plot of the model human compartments.
#'
#' @title Plot human states
#'
#' @param df The dataframe with the data to plot.
#'
#' @export
plot_compartments <- function(df){

  compart_names <- c("Susceptibles", "Infectious", "Recovered")

  time <- max(df$time)

  # 50 year horizon
  brks <- seq(from = 0, to = time, by = 364 * 5)

  ret <- ggplot2::ggplot(df,
                         ggplot2::aes(x = time, y = .data$value, col = .data$diagnostic)) +
    ggplot2::geom_line(size = 0.7) +
    ggplot2::scale_colour_manual(
      values = c("#3333FF", "#FFA500", "#CC0000", "#339900"),
      name = NULL,
      labels = compart_names) +
    ggplot2::scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
    ggplot2::scale_y_continuous(name ="Proportion of population") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = element_text(size = 12),
                   axis.title.y = element_text(size = 12),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12),
                   legend.title = element_text(size = 12),
                   legend.text = element_text(size = 12),
                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

  ret

}


# -----------------------------------------------------------------------------

#' The function makes a faceted plot of the model demographics outputs.
#'
#' @title Plot demographics
#'
#' @param df The dataframe with the data to plot.
#'
#' @export
plot_demographics <- function(df){

  number_of_plots_per_page <- 4

  diagno_nms <- levels(df$diagnostic)

  total_no_plots <- length(diagno_nms)

  no_pages <- ceiling(total_no_plots / number_of_plots_per_page)

  time <- max(df$time)

  brks <- seq(from = 0, to = time, by = 364 * 5)

  # xstrips_labs <- as_labeller(setNames(diagno_nms, levels(df$diagnostics)))

  out <- list()

  for (i in seq_len(no_pages)){

    ret <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = .data$value)) +
      ggplot2::geom_line(size = 0.5, colour = "#63B8FF") +
      ggforce::facet_wrap_paginate(~ diagnostic,
                                   ncol = 2,
                                   nrow = 2,
                                   scales = "free_y",
                                   # labeller = xstrips_labs,
                                   page = i) +
      ggplot2::scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
      ggplot2::scale_y_continuous(name = "", labels = scales::comma) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title.x = element_text(size = 12),
                     axis.title.y = element_text(size = 12),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     strip.text.x = element_text(size = 10),
                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

    out[[i]] <- ret

  }

  out

}



# -----------------------------------------------------------------------------

#' The function makes plots of the mosquito larvae carrying capacity, virus extrinsic
#'  incubation period, adult mosquito mortality rate and reproduction number,
#'  averaged across patches, and arrange them in a grid.
#'
#' @title Plot seasonally-varying variables
#'
#' @param out The dataframe with the data to plot.
#'
#' @importFrom ggplot2 aes geom_line scale_x_continuous
#'   scale_y_continuous ggplot theme_bw
#'
#' @importFrom gridExtra arrangeGrob
#'
#' @export


plot_Kc_eip_delta <- function(out) {

  time <- max(out$TIME)

  my_breaks <- seq(from = 0, to = time, by = 364 * 5)

  Kc_r <- data.frame(x = out$TIME, y = out$Kcav)

  Kc_p <- ggplot(data = Kc_r, aes(x = Kc_r$x, y = Kc_r$y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch Kc") +
    ggtitle("Carrying capacity") +
    theme_bw()

  eip_r <- data.frame(x = out$TIME, y = out$eipav)

  eip_p <- ggplot(data = eip_r, aes(x = eip_r$x, y = eip_r$y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch EIP") +
    ggtitle("Extrinsic Incubation Period") +
    theme_bw()

  delta_r <- data.frame(x = out$TIME, y = out$Deltaav)

  delta_p <- ggplot(data = delta_r, aes(x = delta_r$x, y = delta_r$y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch Delta") +
    ggtitle("Adult mosquito daily mortality rate") +
    theme_bw()

  R0_r <- data.frame(x = out$TIME, y = out$R0t_1av)

  browser()

  R0_p <- ggplot(data = R0_r, aes(x = R0_r$x, y = R0_r$y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous(expression('Mean patch R'['0']*'')) +
    ggtitle(expression('Reproduction number, R'['0']*'')) +
    theme_bw()

  arrangeGrob(Kc_p, eip_p, delta_p, R0_p, ncol = 2)

}
