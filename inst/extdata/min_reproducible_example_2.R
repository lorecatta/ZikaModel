na <- user()
NP <- user()

init_S[,,] <- user() # 3, 2, 4
agert[] <- user()    # 3
deathrt[] <- user()  # 3
rho1[] <- user()     # 2
FOI1[] <- user()     # 4

initial(S[,,]) <- init_S[i,j,k]

# dim(agert_array) <- c(na, 2, NP)
# agert_array[,,] <- 1
#
# dim(deathrt_array) <- c(na, 2, NP)
# deathrt_array[,,] <- 1
#
# dim(rho1_array) <- c(na, 2, NP)
# rho1_array[,,] <- 1
#
# dim(FOI1_array) <- c(na, 2, NP)
# FOI1_array[,,] <- 1
#
# agert_array[1:na,,] <- agert[i]
# deathrt_array[1:na,,] <- deathrt[i]
# rho1_array[,1:2,] <- rho1[j]
# FOI1_array[,,1:NP] <- FOI1[k]
#
# O_S_prob[,,] <- rho1_array[i,j,k] * FOI1_array[i,j,k] + deathrt_array[i,j,k] + agert_array[i,j,k]

O_S_prob[,,] <- rho1[j] * FOI1[k] + agert[i] + deathrt[i]
O_S[,,] <- S[i,j,k] * O_S_prob[i,j,k]
update(S[,,]) <- S[i,j,k] - O_S[i,j,k]

dim(init_S) <- c(na, 2, NP)
dim(rho1) <- 2
dim(FOI1) <- NP
dim(agert) <- na
dim(deathrt) <- na
dim(S) <- c(na, 2, NP)
dim(O_S) <- c(na, 2, NP)
dim(O_S_prob) <- c(na, 2, NP)

output(O_S_prob[,,]) <- TRUE
#output(rho1_array[,,]) <- TRUE
#output(FOI1_array[,,]) <- TRUE
#output(deathrt_array[,,]) <- TRUE
#output(agert_array[,,]) <- TRUE
