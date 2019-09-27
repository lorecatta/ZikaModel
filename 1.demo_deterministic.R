
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

time <- 364 * 150

mpl <- model_param_list_create(season = TRUE)

# generate initial state variables from equilibrium solution
state_init <- equilibrium_init_create(agec = agec,
                                      death = death,
                                      nn_links = nn_links,
                                      model_parameter_list = mpl)

odin_model_path <- system.file("extdata/odin_model_determ.R", package = "ZikaModel")

gen <- odin::odin(odin_model_path,verbose = FALSE)

state_use <- state_init[names(state_init) %in% names(formals(gen))]

mod <- gen(user = state_use)
tt <- seq(0, time, 1)

mod_run <- mod$run(tt)

# shape output
out <- mod$transform_variables(mod_run)

diagno_hum <- c("S", "I1", "R1", "births", "inf_1", "Y1T", "infectious1", "O_S", "D_S", "A_S")

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
df_H$wIR_inf <- lag_diff(df_H$inf_1_cum, 7)

df_H$wIR_inf <- df_H$wIR_inf / df_H$Nt * 1000

df_H$time <- tt
df_H_melt <- melt(df_H,
                  id.vars = "time",
                  variable.name = "diagnostic")

diagno_levs <- c("S", "I1", "R1", "Nt", "births", "inf_1", "Y1T", "infectious1", "O_S", "D_S", "A_S", "inf_1_cum", "Sp", "I1p", "R1p", "wIR_inf")

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

# sum across ages and vaccine status (dims 2 and 3)
births_patch_df <- as.data.frame(dia_hum$births)
names(births_patch_df) <- seq_len(21)
births_patch_df$time <- tt
births_patch_df_melt <- melt(births_patch_df,
                             id.vars = "time",
                             variable.name = "patch")

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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

deaths_patch <- apply(dia_hum$D_S, c(1, 4), sum)
deaths_patch_df <- as.data.frame(deaths_patch)
names(deaths_patch_df) <- seq_len(21)
deaths_patch_df$time <- tt
deaths_patch_df_melt <- melt(deaths_patch_df,
                             id.vars = "time",
                             variable.name = "patch")

brks <- seq(from = 0, to = time, by = 364*10)
p <- ggplot(deaths_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "deaths (S)") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "deaths_S_by_patch.png",
          wdt = 15,
          hgt = 15)

aged_patch <- apply(dia_hum$A_S, c(1, 4), sum)
aged_patch_df <- as.data.frame(aged_patch)
names(aged_patch_df) <- seq_len(21)
aged_patch_df$time <- tt
aged_patch_df_melt <- melt(aged_patch_df,
                             id.vars = "time",
                             variable.name = "patch")

brks <- seq(from = 0, to = time, by = 364*10)
p <- ggplot(aged_patch_df_melt) +
  geom_line(aes(x = time, y = value), colour = "#63B8FF") +
  ggplot2::facet_wrap(~ patch, ncol = 4) +
  scale_y_continuous(name = "aged (S)") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir_1,
          out_fl_nm = "aged_S_by_patch.png",
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

brks <- seq(from = 0, to = time, by = 364*10)
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

brks <- seq(from = 0, to = time, by = 364*10)
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
O_S_prob_full_melt$age <- factor(O_S_prob_full_melt$age,
                                 levels = unique(O_S_prob_full_melt$age),
                                 labels = unique(O_S_prob_full_melt$age))

O_S_prob_melt <- subset(O_S_prob_full_melt, vaccine == 1 & patch == 1)

brks <- seq(from = 0, to = time, by = 364*10)
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
S_full_melt$age <- factor(S_full_melt$age,
                          levels = unique(S_full_melt$age),
                          labels = unique(S_full_melt$age))

S_melt <- subset(S_full_melt, vaccine == 1 & patch == 1)

brks <- seq(from = 0, to = time, by = 364*10)
p <- ggplot(S_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "S") +
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

D_S_full_melt <- melt(out$D_S)
names(D_S_full_melt) <- c("time", "age", "vaccine", "patch", "value")
D_S_full_melt$age <- factor(D_S_full_melt$age,
                            levels = unique(D_S_full_melt$age),
                            labels = unique(D_S_full_melt$age))

D_S_melt <- subset(D_S_full_melt, vaccine == 1 & patch == 1)

brks <- seq(from = 0, to = time, by = 364*10)
p <- ggplot(D_S_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "D_S") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("D_S_vaccine_%s_patch_%s.png", 1, 1),
          wdt = 10,
          hgt = 8)

A_S_full_melt <- melt(out$A_S)
names(A_S_full_melt) <- c("time", "age", "vaccine", "patch", "value")
A_S_full_melt$age <- factor(A_S_full_melt$age,
                            levels = unique(A_S_full_melt$age),
                            labels = unique(A_S_full_melt$age))

A_S_melt <- subset(A_S_full_melt, vaccine == 1 & patch == 2)

brks <- seq(from = 0, to = time, by = 364*10)
p <- ggplot(A_S_melt) +
  geom_line(aes(x = time, y = value, colour = age)) +
  scale_fill_viridis() +
  scale_y_continuous(name = "A_S") +
  scale_x_continuous(name = "Years", breaks = brks, labels = brks / 364) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8))
save_plot(p,
          out_dir,
          out_fl_nm = sprintf("A_S_vaccine_%s_patch_%s.png", 1, 2),
          wdt = 10,
          hgt = 8)


# model_run <- run_model(agec,
#                        death,
#                        nn_links,
#                        time = time,
#                        season = TRUE)
#
# save_plot(model_run$plot,
#           out_dir,
#           "compartments_human",
#           wdt = 17,
#           hgt = 12)
#
# post_processing(model_run$dat, out_dir)
