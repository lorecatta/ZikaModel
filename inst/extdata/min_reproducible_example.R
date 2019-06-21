na <- user()
NP <- user()

dim(S) <- c(na, 2, NP)
dim(init_S) <- c(na, 2, NP)
init_S[,,] <- user()
initial(S[,,]) <- init_S[i,j,k]

dim(age_S_prob) <- na
age_S_prob[] <- user()
dim(age_S) <- c(na, 2, NP)

age_S[,,] <- S[i,j,k] * age_S_prob[i]

update(S[,,]) <- S[i,j,k] - age_S[i,j,k]

dim(ST_sum) <- c(na,NP)
ST_sum[1,1:NP] <- S[i,1,j] + S[i,2,j]
ST_sum[2:na,1:NP] <- ST_sum[i-1,j] + S[i,1,j] + S[i,2,j]

#dim(ST) <- NP
ST <- ST_sum[na,1] / 10

incub <- user()
inf_per <- user()
dim(infectious1) <- NP
initial(infectious1[]) <- 0

update(infectious1[]) <- Y1T_del_inc + infectious1[i] - Y1T_del_inc_ip

Y1T_del_inc <- delay(ST, 0)
#dim(Y1T_del_inc) <- NP
del_amount <- inf_per + incub
Y1T_del_inc_ip <- delay(ST, del_amount)
#dim(Y1T_del_inc_ip) <- NP

dim(test) <- 10
initial(test[]) <- runif(-1, 2)
update(test[]) <- max(min(test[i], 1), 0)

output(ST) <- TRUE
output(Y1T_del_inc) <- TRUE
output(Y1T_del_inc_ip) <- TRUE
