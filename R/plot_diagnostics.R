#------------------------------------------------
#' plot_diagnostics
#'
#' \code{plot_diagnostics} makes a faceted plot of the model diagnostics.
#'
#' @param df The dataframe with the data to plot.
#' @param diagno_nms A character string of descriptive names of the diagnostics. This will be used for the facet y labels.
#'
#' @inheritParams save_plot
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   scale_y_continuous xlab ylab geom_line ggplot theme_bw
#'   theme element_text margin unit ggtitle as_labeller
#'
#' @export


plot_diagnostics <- function(df, out_pth, out_fl_nm, diagno_nms){

  number_of_plots_per_page <- 4

  total_no_plots <- length(diagno_nms)

  no_pages <- ceiling(total_no_plots / number_of_plots_per_page)

  time <- max(df$time)

  brks <- seq(from = 0, to = time, by = 364 * 5)

  # xstrips_labs <- as_labeller(setNames(diagno_nms, levels(df$diagnostics)))

  dir.create(out_pth, FALSE, TRUE)

  plot_ttl <- gsub(".*_", "", out_fl_nm)

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
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            plot.title = element_text(margin = margin(0,0,0.5,0,"cm"))) +
      ggtitle(paste0("SEIR Zika model - ", plot_ttl))

    save_plot(ret,
              out_pth,
              out_fl_nm = sprintf("%s_%s%s", out_fl_nm, i, ".png"),
              wdt = 18,
              hgt = 10)

  }

}
