
devtools::load_all()

out_dir <- file.path("figures", "deterministic_no_seasonality")

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

integer_time_steps <- 364 * 100

mpl <- model_param_list_create(season = FALSE)

# generate initial state variables from equilibrium solution
state_init <- equilibrium_init_create(agec = agec,
                                      death = death,
                                      nn_links = nn_links,
                                      model_parameter_list = mpl)

odin_model_path <- system.file("extdata/odin_model_determ.R", package = "ZikaModel")

gen <- odin::odin(odin_model_path,verbose = FALSE)

state_use <- state_init[names(state_init) %in% names(formals(gen))]

mod <- gen(user = state_use)

its <- seq(0, integer_time_steps, 1)

mod_run <- mod$run(its)

# shape output
out <- mod$transform_variables(mod_run)

tt <- out$TIME
time <- max(tt)

diagno_hum <- c("S", "I1", "R1", "births", "inf_1", "Y1T", "infectious1", "O_S")

dia_hum <- setNames(out[diagno_hum], diagno_hum)

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

mat_H <- cbind(mat_H, prop)

df_H <- as.data.frame(mat_H)

# rate of total weekly infections
df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, 14)

df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

df_H$time <- tt
df_H_melt <- melt(df_H,
                  id.vars = "time",
                  variable.name = "diagnostic")

diagno_levs <- c("S", "I1", "R1", "Nt", "births", "inf_1", "Y1T", "infectious1", "O_S", "inf_1_cum", "Sp", "I1p", "R1p", "wIR_inf")

df_H_melt$diagnostic <- factor(df_H_melt$diagnostic, levels = diagno_levs, labels = diagno_levs)

plot_diagnostics(df_H_melt,
                 out_dir,
                 "human_diagnostics",
                 diagno_levs)


# -----------------------------------------------------------------------------
#
# Plot diagnostics by patch
#
# -----------------------------------------------------------------------------


out_dir_1 <- file.path(out_dir, "patch")
brks <- seq(from = 0, to = time, by = 364 * 5)

# sum across ages and vaccine status (dims 2 and 3)
births_patch_df <- as.data.frame(dia_hum$births)
names(births_patch_df) <- seq_len(21)
births_patch_df$time <- tt
births_patch_df_melt <- melt(births_patch_df,
                             id.vars = "time",
                             variable.name = "patch")

p <- ggplot(births_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "births") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "births_by_patch.png",
          wdt = 15,
          hgt = 15)

# sum across ages and vaccine status (dims 2 and 3)
inf_1_patch <- apply(dia_hum$inf_1, c(1, 4), sum)
inf_1_patch_df <- as.data.frame(inf_1_patch)
names(inf_1_patch_df) <- seq_len(21)
inf_1_patch_df$time <- tt
inf_1_patch_df_melt <- melt(inf_1_patch_df,
                         id.vars = "time",
                         variable.name = "patch")

p <- ggplot(inf_1_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "inf_1") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "inf_1_by_patch.png",
          wdt = 15,
          hgt = 15)

# sum across ages and vaccine status (dims 2 and 3)
S_patch <- apply(dia_hum$S, c(1, 4), sum)
S_patch_df <- as.data.frame(S_patch)
names(S_patch_df) <- seq_len(21)
S_patch_df$time <- tt
S_patch_df_melt <- melt(S_patch_df,
                        id.vars = "time",
                        variable.name = "patch")

p <- ggplot(S_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "S") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "S_by_patch.png",
          wdt = 15,
          hgt = 15)

I1_patch <- apply(dia_hum$I1, c(1, 4), sum)
I1_patch_df <- as.data.frame(I1_patch)
names(I1_patch_df) <- seq_len(21)
I1_patch_df$time <- tt
I1_patch_df_melt <- melt(I1_patch_df,
                        id.vars = "time",
                        variable.name = "patch")

p <- ggplot(I1_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "I1") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "I1_by_patch.png",
          wdt = 15,
          hgt = 15)

R1_patch <- apply(dia_hum$R1, c(1, 4), sum)
R1_patch_df <- as.data.frame(R1_patch)
names(R1_patch_df) <- seq_len(21)
R1_patch_df$time <- tt
R1_patch_df_melt <- melt(R1_patch_df,
                         id.vars = "time",
                         variable.name = "patch")

p <- ggplot(R1_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "R1") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "R1_by_patch.png",
          wdt = 15,
          hgt = 15)

FOI1p_patch_df <- as.data.frame(out$FOI1p)
names(FOI1p_patch_df) <- seq_len(21)
FOI1p_patch_df$time <- tt
FOI1p_patch_df_melt <- melt(FOI1p_patch_df,
                         id.vars = "time",
                         variable.name = "patch")

p <- ggplot(FOI1p_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "FOI1p") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "FOI1p_by_patch.png",
          wdt = 15,
          hgt = 15)


# -----------------------------------------------------------------------------
#
# Plot diagnostics by age group
#
# -----------------------------------------------------------------------------


out_dir_2 <- file.path(out_dir, "age")

deathrt_age_df <- as.data.frame(out$deathrt)
names(deathrt_age_df) <- seq_len(11)
deathrt_age_df$time <- tt
deathrt_age_df_melt <- melt(deathrt_age_df,
                             id.vars = "time",
                             variable.name = "age")

p <- ggplot(deathrt_age_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ age, ncol = 4#,
                      # scales = "free_y"
                      ) +
  scale_y_continuous(name = "deathrt") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8)) +
  ggplot2::coord_cartesian(ylim=c(0,.004))
save_plot(p,
          out_dir_2,
          out_fl_nm = "deathrt_by_age.png",
          wdt = 15,
          hgt = 15)

agert_age_df <- as.data.frame(out$agert)
names(agert_age_df) <- seq_len(11)
agert_age_df$time <- tt
agert_age_df_melt <- melt(agert_age_df,
                            id.vars = "time",
                            variable.name = "age")

p <- ggplot(agert_age_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ age, ncol = 4, scales = "free_y") +
  scale_y_continuous(name = "agert") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_2,
          out_fl_nm = "agert_by_age.png",
          wdt = 15,
          hgt = 15)



# -----------------------------------------------------------------------------
#
# Plot diagnostics by all
#
# -----------------------------------------------------------------------------



library(viridis)

O_S_prob_full_melt <- melt(out$O_S_prob)
names(O_S_prob_full_melt) <- c("time", "age", "vaccine", "patch", "value")
no_age <- length(unique(O_S_prob_full_melt$age))
no_vaccine <- length(unique(O_S_prob_full_melt$vaccine))
no_patch <- length(unique(O_S_prob_full_melt$patch))
combs <- no_age * no_vaccine * no_patch
tt_long <- rep(tt, combs)
O_S_prob_full_melt$time <- tt_long
O_S_prob_full_melt$age <- factor(O_S_prob_full_melt$age,
                                 levels = unique(O_S_prob_full_melt$age),
                                 labels = unique(O_S_prob_full_melt$age))

O_S_prob_melt <- subset(O_S_prob_full_melt, vaccine == 1 & patch == 1)

p <- ggplot(O_S_prob_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "O_S_prob") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("O_S_prob_vaccine_%s_patch_%s.png", 1, 1),
          wdt = 10,
          hgt = 8)


S_full_melt <- melt(out$S)
names(S_full_melt) <- c("time", "age", "vaccine", "patch", "value")
S_full_melt$time <- tt_long
S_full_melt$age <- factor(S_full_melt$age,
                          levels = unique(S_full_melt$age),
                          labels = unique(S_full_melt$age))

S_melt <- subset(S_full_melt, vaccine == 1 & patch == 1)

p <- ggplot(S_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "S", breaks = seq(0,6e+6,1e+6), labels = seq(0,6e+6,1e+6), limits = c(0,6e+6)) +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("S_vaccine_%s_patch_%s.png", 1, 1),
          wdt = 10,
          hgt = 8)

I1_full_melt <- melt(out$I1)
names(I1_full_melt) <- c("time", "age", "vaccine", "patch", "value")
I1_full_melt$time <- tt_long
I1_full_melt$age <- factor(I1_full_melt$age,
                           levels = unique(I1_full_melt$age),
                           labels = unique(I1_full_melt$age))

I1_melt <- subset(I1_full_melt, vaccine == 1 & patch == 1)

p <- ggplot(I1_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "I1", breaks = seq(0,1.4e+6,2e+5), labels = seq(0,1.4e+6,2e+5), limits = c(0,1.4e+6)) +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("I1_vaccine_%s_patch_%s.png", 1, 1),
          wdt = 10,
          hgt = 8)

O_S_full_melt <- melt(out$O_S)
names(O_S_full_melt) <- c("time", "age", "vaccine", "patch", "value")
O_S_full_melt$time <- tt_long
O_S_full_melt$age <- factor(O_S_full_melt$age,
                           levels = unique(O_S_full_melt$age),
                           labels = unique(O_S_full_melt$age))

O_S_melt <- subset(O_S_full_melt, vaccine == 1 & patch == 1)

p <- ggplot(O_S_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "O_S") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("O_S_vaccine_%s_patch_%s.png", 1, 1),
          wdt = 10,
          hgt = 8)

inf_1_full_melt <- melt(out$inf_1)
names(inf_1_full_melt) <- c("time", "age", "vaccine", "patch", "value")
inf_1_full_melt$time <- tt_long
inf_1_full_melt$age <- factor(inf_1_full_melt$age,
                            levels = unique(inf_1_full_melt$age),
                            labels = unique(inf_1_full_melt$age))

inf_1_melt <- subset(inf_1_full_melt, vaccine == 1 & patch == 3)

p <- ggplot(inf_1_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "inf_1", breaks = seq(0,7.5e+4,1.5e+4), labels = seq(0,7.5e+4,1.5e+4), limits = c(0,7.5e+4)) +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("inf_1_vaccine_%s_patch_%s.png", 1, 3),
          wdt = 10,
          hgt = 8)

model_run <- ZikaModel::run_model(agec,
                                  death,
                                  nn_links,
                                  time = time,
                                  season = TRUE)
#
# save_plot(model_run$plot,
#           out_dir,
#           "compartments_human",
#           wdt = 17,
#           hgt = 12)
#
# post_processing(model_run$dat, out_dir)
