#------------------------------------------------
#' plot_compartments
#'
#' \code{plot_compartments} makes a plot of the model compartments
#'
#' @importFrom ggplot2 aes scale_colour_manual scale_x_continuous
#'   scale_y_continuous xlab ylab geom_line ggplot theme_bw
#'   theme element_text margin unit ggtitle
#'
#' @export


plot_compartments <- function(df, compart_names, ttl, out_pth = NULL, out_fl_nm = NULL){

  brks <- seq(from = 0, to = time, by = 364*5)

  ret <- ggplot(df,
                aes(x = time, y = .data$value, col = .data$compartment)) +
    geom_line(size = 0.7) +
    # facet_wrap(~ .data$vaccine, nrow = 2, labeller = xstrips_labs) +
    scale_colour_manual(
      values = c("#3333FF", "#FFA500", "#CC0000", "#339900"),
      name = NULL,
      labels = compart_names) +
    scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
    scale_y_continuous(name ="Proportion of population") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.title = element_text(margin = margin(0,0,0.5,0,"cm"))) +
    ggtitle(ttl)

  if(!is.null(out_pth)){

    dir.create(out_pth, FALSE, TRUE)
    png(file.path(out_pth, out_fl_nm),
        width = 17,
        height = 12,
        units = "cm",
        pointsize = 12,
        res = 300)
    print(ret)
    dev.off()

  } else {

    ret

  }

}
