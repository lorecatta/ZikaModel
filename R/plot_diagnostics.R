#------------------------------------------------
#' plot_demographics
#'
#' \code{plot_demographics} makes a faceted plot of the model diagnostics.
#'
#' @param df The dataframe with the data to plot.
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   scale_y_continuous xlab ylab geom_line ggplot theme_bw
#'   theme element_text margin unit ggtitle as_labeller
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

    ret <- ggplot(df, aes(x = time, y = value)) +
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
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            strip.text.x = element_text(size = 8),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

    out[[i]] <- ret

  }

  out

}
