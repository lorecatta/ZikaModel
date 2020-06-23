## code to prepare `DATASET` dataset goes here


# patch structure and seasonality ---------------------------------------------


nn_links <- as.matrix(read.table(file.path("data-raw", "nn_links.txt"),
                                 check.names = FALSE,
                                 header = TRUE,
                                 row.names = 1))

usethis::use_data(nn_links)

amplitudes_phases <- as.matrix(read.csv(file.path("data-raw", "amplitudes_phases.csv")))

usethis::use_data(amplitudes_phases)


# microcephaly cases ----------------------------------------------------------


# from Hay et al. 2018. Potential inconsistencies in Zika surveillance data and
# our understanding of risk during pregnancy. PLOSNTD

# Reads in all files ending "multivariate_chain" from specified working dir
chains <- zikaInfer::load_mcmc_chains(location = "data-raw/hay_et_al_2018",
                                      asList = FALSE,
                                      unfixed = FALSE)

## Gamma function used
# zikaInfer::microceph_v1

chain_means <- colMeans(chains)
gamma_mean <- as.numeric(chain_means["mean"])
gamma_var <- as.numeric(chain_means["var"])
gamma_c <- as.numeric(chain_means["c"])

pars <- c(mean = gamma_mean,
          var = gamma_var,
          c = gamma_c)

# probability that a foetus developed microcephaly
# given that the mother was infected in a particular week during pregnancy
m_risk_probs <- zikaInfer::microceph_v1(pars = pars)

mc_prob_ZIKV_pregn <- data.frame(day = 0:279, prob = m_risk_probs)

usethis::use_data(mc_prob_ZIKV_pregn)


# birth rates for Brazil ------------------------------------------------------


# from World Population Prospects 2019
# (https://population.un.org/wpp/Download/Standard/Fertility/)
# birth rates for Brazil for age classes:
# 15-19, 20-24, 25-29, 30-34, 35-39, 40-44 and 45-49.

na <- 11

br_brazil_raw <- read.csv(file.path("data-raw", "age_specific_birth_rates_Brazil.csv"))

br_brazil <- br_brazil_raw[, "X2010...2015"] / 1000

br_brazil_2 <- c(0, br_brazil)

br_brazil_av <- c()

my_seq <- seq(1, length(br_brazil_2), by = 2)

for (i in seq_along(my_seq)) {

  index <- my_seq[i]

  out <- (br_brazil_2[index] + br_brazil_2[index + 1]) / 2

  br_brazil_av[i] <- out

}

br_brazil_age <- rep(0, na)

br_brazil_age[1] <- 0
br_brazil_age[2] <- 0
br_brazil_age[3] <- br_brazil_av[1]
br_brazil_age[4] <- br_brazil_av[2]
br_brazil_age[5] <- br_brazil_av[3]
br_brazil_age[6] <- br_brazil_av[4]
br_brazil_age[7] <- 0
br_brazil_age[8] <- 0
br_brazil_age[9] <- 0
br_brazil_age[10] <- 0
br_brazil_age[11] <- 0

# br_brazil_age_2 <- br_brazil_age / YL * DT

usethis::use_data(br_brazil_age)
