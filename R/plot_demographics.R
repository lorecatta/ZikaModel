
# -----------------------------------------------------------------------------

#' The function makes a faceted plot of the model demographics outputs.
#'
#' @title Plot demographics
#'
#' @param df The dataframe with the data to plot.
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   scale_y_continuous xlab ylab geom_line ggplot theme_bw
#'   theme element_text margin unit ggtitle .data as_labeller
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

    ret <- ggplot(df, aes(x = time, y = .data$value)) +
      geom_line(size = 0.5, colour = "#63B8FF") +
      ggforce::facet_wrap_paginate(~ diagnostic,
                                   ncol = 2,
                                   nrow = 2,
                                   scales = "free_y",
                                   # labeller = xstrips_labs,
                                   page = i) +
      scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
      scale_y_continuous(name = "", labels = scales::comma) +
      theme_bw() +
      theme(axis.title.x = element_text(size = 12),
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
#'  incubation period and adult mosquito mortality rate, and arrange them in a grid.
#'  These variables are those directly affected by seasonality in the model.
#'
#' @title Plot seasonally-varying variables
#'
#' @param df The dataframe with the data to plot.
#'
#' @importFrom ggplot2 aes geom_line scale_x_continuous
#'   scale_y_continuous ggplot theme_bw
#'
#' @export


plot_Kc_eip_delta <- function(out) {

  time <- max(out$TIME)

  my_breaks <- seq(from = 0, to = time, by = 364 * 5)

  Kc_r <- data.frame(x = out$TIME, y = out$Kcav)

  Kc_p <- ggplot(data = Kc_r, aes(x = x, y = y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch Kc") +
    ggtitle("Carrying capacity") +
    theme_bw()

  eip_r <- data.frame(x = out$TIME, y = out$eipav)

  eip_p <- ggplot(data = eip_r, aes(x = x, y = y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch EIP") +
    ggtitle("Extrinsic Incubation Period") +
    theme_bw()

  delta_r <- data.frame(x = out$TIME, y = out$Deltaav)

  delta_p <- ggplot(data = delta_r, aes(x = x, y = y)) +
    geom_line(color = 'royalblue', size = 0.5) +
    scale_x_continuous("Years", breaks = my_breaks, labels = my_breaks/364) +
    scale_y_continuous("Mean patch Delta") +
    ggtitle("Adult mosquito daily mortality rate") +
    theme_bw()

  arrangeGrob(Kc_p, eip_p, delta_p, ncol = 2)

}
