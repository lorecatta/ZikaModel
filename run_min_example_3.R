gen <- odin::odin("inst/extdata/min_reproducible_example_3.R", verbose = FALSE)

na <- 10
init_S <- seq_len(na)

pars <- list(na = na,
             init_S = init_S)

mod <- gen(user = pars)

t <- seq_len(3)

y <- mod$run(t)

out <- mod$transform_variables(y)
