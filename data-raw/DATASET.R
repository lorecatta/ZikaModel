## code to prepare `DATASET` dataset goes here

nn_links <- read.table(file.path("data-raw", "nn_links.txt"), check.names = FALSE, header = TRUE)

usethis::use_data(nn_links)
