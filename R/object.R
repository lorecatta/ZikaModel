
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
#' @param keep name of variable to stratify by
#'   (allowed are \code{c("patch", "vaccine")} for type H and
#'   \code{"patch"} for type M. Default is no stratification)
#'
#' @param ... additional arguments affecting the plot produced.
#'
#' @export
plot.Zika_model_simulation <- function(x,
                                       type,
                                       var_select = NULL,
                                       keep = NULL,
                                       ...) {

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
      ggplot2::facet_wrap(stats::as.formula(paste("~", keep)),
                          ncol = 4,
                          scales = "free_y")
  } else if (keep == "all") {
    pds <- subset(pds, patch == 1)
    p <- p + ggplot2::geom_line(data = pds,
                                ggplot2::aes(x = .data$t, y = .data$y,
                                             col = .data$age)) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", "vaccine")),
                          ncol = 1,
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
