
gen <- odin::odin("inst/extdata/min_reproducible_example_2.R", verbose = FALSE)

FOI1 <- c(0.01,
          0.02,
          0.03,
          0.04)

agert <- c(1, 2, 3)

deathrt <- c(10, 20, 30)

rho1 <- c(1, 0)

pars <- list(na = 3,
             NP = 4,
             init_S = array(rep(1, 24), c(3, 2, 4)),
             rho1 = rho1,
             FOI1 = FOI1,
             agert = agert,
             deathrt = deathrt)

mod <- gen(user = pars)

t <- seq_len(3)

y <- mod$run(t)

out <- mod$transform_variables(y)
