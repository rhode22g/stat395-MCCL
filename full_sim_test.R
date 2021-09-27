# date created: September 1, 2021

# This script runs the process of simulation for the full ancestor dist'n estimation for seed = 35 as a test value

######
# loading needed packages
library(ggplot2)
library(sfsmisc)
library(xlsx)
library(writexl)
######

plot(x = ancestor_sim1[,1], type = "h")

######
# FIND TRUE VALUES
an_seq <- matrix(c(0:(2^(20)-1)), ncol=1)
an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=20) ) 

#true values for sites 1-19 phi_s
true_phi_sim1 <- rep(NA, 19)
for (j in 1:19){
  to_sum <- rep(NA, 2^20)
  for (h in 1:(2^20)){
    if (an_seq[h, j] == 0 & an_seq[h, (j+1)]==0){
      to_sum[h] <- ancestor_sim1[h]
    }
  }
  true_phi_sim1[j] <- sum(to_sum, na.rm = TRUE)
}
true_phi_sim1

# number of ancestor sequences that are non-zero
num_zero <- sum(ancestor_sim1[,1] == 0)
num_nonzero <- 2^(20) - num_zero
num_nonzero
######

######
# initialize matrix to store bias and SD
test_simulation_35 <- matrix(NA, nrow = 1, ncol = 19)
colnames(test_simulation_35) <- paste("site_",seq(1:19),sep="")
test_full_35 <- matrix(NA, nrow=num_nonzero, ncol= 3)
colnames(test_full_35) <- c("Bitwise", "True", "Estimate")
indx <- 1
for (i in 1:dim(ancestor_sim1)[1]){
  if (ancestor_sim1[i,1] != 0){
    test_full_35[indx,1] <- i
    test_full_35[indx,2] <- ancestor_sim1[i,1]
    indx <- indx + 1
  }
}
# test_full_35[,1] <- ancestor_sim1[,1]*(ancestor_sim1[,1] != 0)
# test_full_35[,1:2] <- t( apply(ancestor_sim1, ((ancestor_sim1[,1] != 0))))
test_full_35 <- as.data.frame(test_full_35)
write_xlsx(test_full_35, path = "C:/Users/gbean/Documents/MCCL_full_testseed_m1")
# write.xlsx(test_full_35, file = "C:/Users/gbean/Documents/MCCL_full_testseed_m1","Test Seed", append = TRUE)

# ONE REPLICATION AS TEST CASE
# find descendant sample using seed = 35
descend_sample_35 <- descendent_sample(L=20, q=0.005, n=100, seed = 35, hapmap = hapmap_binary)
  
# find onewise estimates
onewise_estimates_35 <- estimates_m0(descendents = descend_sample_35, L=20, n=100)
  
# find pairwise estimates
pairwise_estimates_35 <- estimates_m1(L=20, q=0.005, n=100, d= descend_sample_35)

# find full estimates
full_estimates_35 <- ancestor_pairs(L=20, pairs = pairwise_estimates_35, ones = onewise_estimates_35)
  
# append estimate for this replication into the matrix
test_simulation_35[1,] <- pairwise_estimates_35[1,]
test_full_35[,3] <- full_estimates_35[(test_full_35[,1]),1]

# calculate performance measures
prop_correct_nonzeros <- (sum(test_full_35[,3] != 0 & test_full_35[,2] != 0)) / num_nonzero
prop_correct_nonzeros

amnt_incorrect_prob <- 1 - sum(test_full_35[,3])
amnt_incorrect_prob

bias_35 <- vector(length = num_nonzero)
bias_35 <- test_full_35[,3] - test_full_35[,2]
bias_35
plot(x = 1:num_nonzero , y = bias_35)
#####

####
# send to excel
test_full_35 <- as.data.frame(test_full_35)
write_xlsx(test_full_35, path = "C:/Users/gbean/Documents/MCCL_full_testseed_m1.xlsx")
####
