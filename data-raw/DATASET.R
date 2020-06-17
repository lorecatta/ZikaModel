## code to prepare `DATASET` dataset goes here

nn_links <- as.matrix(read.table(file.path("data-raw", "nn_links.txt"),
                                 check.names = FALSE,
                                 header = TRUE,
                                 row.names = 1))

amplitudes_phases <- as.matrix(read.csv(file.path("data-raw", "amplitudes_phases.csv")))

usethis::use_data(nn_links)
usethis::use_data(amplitudes_phases)

# from James Hay
chain_1 <- read.csv(file.path("data-raw",
                              "hay_et_al_2018",
                              "test_1_multivariate_chain.csv"))

chain_2 <- read.csv(file.path("data-raw",
                              "hay_et_al_2018",
                              "test_2_multivariate_chain.csv"))

chain_3 <- read.csv(file.path("data-raw",
                              "hay_et_al_2018",
                              "test_3_multivariate_chain.csv"))

mc_prob_ZIKV_pregn <- calculate_mc_probability_given_ZIKV_infection_pregnancy()

usethis::use_data(mc_prob_ZIKV_pregn)

# from World Population Prospects 2019
# (https://population.un.org/wpp/Download/Standard/Fertility/)

br_brazil_raw <- read.csv(file.path("data", "age_specific_birth_rates_Brazil.csv"))

br_brazil <- calculate_brazil_birth_rates(br_brazil_raw)

usethis::use_data(br_brazil)
