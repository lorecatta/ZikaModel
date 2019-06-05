
devtools::load_all()

out_dir <- file.path("figures", "deterministic")

agec <- c(1, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10)

death <- c(1e-10,
           1e-10,
           1e-10,
           0.00277068683332695,
           0.0210680857689784,
           0.026724997685722,
           0.0525354529367476,
           0.0668013582441452,
           0.119271483740379,
           0.279105747097929,
           0.390197266957464)

time <- 21840 # 18200 # 50 years

odin_model_path <- system.file("extdata/odin_model_determ.R", package = "ZikaModel")

wh <- create_r_model(odin_model_path = odin_model_path,
                     agec = agec,
                     death = death,
                     nn_links = nn_links,
                     season = TRUE)

mod <- wh$generator(user = wh$state)

tt <- seq(0, time, 1)

mod_run <- mod$run(tt)

out <- mod$transform_variables(mod_run)

# shape output
out <- mod$transform_variables(mod_run)
hum <- list("S" = out$S, "I" = out$I1, "R" = out$R1, "Nt" = out$Ntotal)
humsum <- lapply(hum, function(x){apply(x, sum, MARGIN = 1)})
mat_H <- do.call("cbind", humsum[c("S", "I", "R")])
mat_H <- mat_H / humsum$Nt
mat_H[is.na(mat_H)] <- 0
df_H <- as.data.frame(mat_H)
df_H$time <- tt
df_H_melt <- melt(df_H,
                  id.vars = "time",
                  variable.name = "compartment")

ret <- plot_compartments(df_H_melt,
                         c("Susceptibles", "Infectious", "Recovered"),
                         "SEIR Zika model - human states")

save_plot(ret,
          out_dir,
          "compartments_human",
          wdt = 17,
          hgt = 12)

post_processing(out, out_dir)
