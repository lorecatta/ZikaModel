
#------------------------------------------------------------------------------
#' save_plot
#'
#' \code{save_plot} save a png file of a plot
#'
#' @param plot_obj The plot object
#' @param out_pth The path where to save the plot
#' @param out_fl_nm The output file name
#' @param wdt The plot width
#' @param hgt The plot height
#'
#' @export


save_plot <- function(plot_obj, out_pth, out_fl_nm, wdt, hgt){

  dir.create(out_pth, FALSE, TRUE)
  png(file.path(out_pth, paste0(out_fl_nm, ".png")),
      width = wdt,
      height = hgt,
      units = "cm",
      pointsize = 12,
      res = 300)
  print(plot_obj)
  on.exit(dev.off())

}
