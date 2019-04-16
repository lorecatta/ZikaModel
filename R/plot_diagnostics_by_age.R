add_diagno_name_var <- function(x, b){
  b[[x]]$diagnostics <- names(b[x])
  b[[x]]
}

add_time_var <- function(x){
  x$time <- tt
  x
}

give_col_names <- function(x){
  colnames(x) <- as.character(seq_len(11))
  x
}


plot_diagnostics_by_age <- function(df, out_pth, out_fl_nm, diagno_nms){
  
  # browser()
  
  no_pages <- 1
  
  brks <- seq(from = 0, to = time, by = 364*10)
  
  xstrips_labs <- as_labeller(setNames(diagno_nms, levels(df$diagnostics)))
  
  dir.create(out_pth, FALSE, TRUE)
  
  for (i in seq_len(no_pages)){
    
    png(file.path(out_pth, sprintf("%s_%s%s", out_fl_nm, i, ".png")),
        width = 20,
        height = 12,
        units = "cm",
        pointsize = 12,
        res = 300)
    
    ret <- ggplot(df, aes(x = .data$time, y = .data$value, colour = .data$age_group)) +
      geom_line(size = 0.4) +
      ggforce::facet_wrap_paginate(~ .data$diagnostics, 
                                   ncol = 1, 
                                   nrow = 2, 
                                   scales = "free_y", 
                                   labeller = xstrips_labs, 
                                   page = i) +
      scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
      scale_y_continuous(name = "", labels = scales::comma_format(accuracy = .01)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            strip.text.x = element_text(size = 8),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            plot.title = element_text(margin = margin(0,0,0.5,0,"cm"))) + 
      ggtitle("SEIR Zika model - diagnostics")
    
    print(ret)
    
    dev.off()
    
  }
  
}
