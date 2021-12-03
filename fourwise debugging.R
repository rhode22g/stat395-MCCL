sim9_fullb <- matrix(NA, nrow = (2^20), ncol = 102)
sim9_fullb[,1] <- c(1:(2^20))
sim9_fullb[,2] <- ancestor_sim1

# comp_time9 <- proc.time()
# 100 REPLICATIONS
# for (i in 1:100){
  # find descendant sample using index as the seed
i <- 3
descend_sample <- descendent_sample(L=20, q=0.005, n=100, seed = i, hapmap = hapmap_binary)
  
  # find threewise estimates
threewise_estimates <- estimates_m2(L=20, q=0.005, n=100, descend = descend_sample)
  
  # find fourwise estimates
fourwise_estimates <- estimates_m3(L=20, q=0.005, n=100, d = descend_sample)
  
  # append estimate for this replication into the matrix
  # simulation_results9_phi00[i,] <- fourwise_estimates[1,]
  # simulation_results9_phi01[i,] <- fourwise_estimates[3,]
  # simulation_results9_phi10[i,] <- fourwise_estimates[5,]
  # simulation_results9_phi11[i,] <- fourwise_estimates[7,]
  
sim9_fullb[,(i+2)] <- ancestor_four_estimation(L = 20, m = 3, four_est = fourwise_estimates, three_est = threewise_estimates)
# }
  
  
an_seq <- matrix(c(0:(2^(8)-1)), ncol=1)
an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=8) ) 

  
## repeat simulation 9 with another sequence length
## FIND TRUE ANCESTOR
ancestor_sim8 <- ancestor(L=8, hapmap_data = hapmap_binary)


## RUN SIMULATION LOOP
simulation_results9new_phi00 <- matrix(NA, nrow = 100, ncol = 5)
colnames(simulation_results9new_phi00) <- paste("site_",seq(1:5),sep="")
simulation_results9new_phi01 <- matrix(NA, nrow = 100, ncol = 5)
colnames(simulation_results9new_phi01) <- paste("site_",seq(1:5),sep="")
simulation_results9new_phi10 <- matrix(NA, nrow = 100, ncol = 5)
colnames(simulation_results9new_phi10) <- paste("site_",seq(1:5),sep="")
simulation_results9new_phi11 <- matrix(NA, nrow = 100, ncol = 5)
colnames(simulation_results9new_phi11) <- paste("site_",seq(1:5),sep="")

# initialize matrix to store joint estimates on all possible sequences
sim9new_full <- matrix(NA, nrow = (2^8), ncol = 102)
sim9new_full[,1] <- c(1:(2^8))
sim9new_full[,2] <- ancestor_sim8
for (i in 1:100){
  # find descendant sample using index as the seed
  descend_sample9 <- descendent_sample(L=8, q=0.005, n=100, seed = i, hapmap = hapmap_binary)
  
  # find threewise estimates
  threewise_estimates9 <- estimates_m2(L=8, q=0.005, n=100, descend = descend_sample9)
  
  # find fourwise estimates
  fourwise_estimates9 <- estimates_m3(L=8, q=0.005, n=100, d = descend_sample9)
  
  # append estimate for this replication into the matrix
  simulation_results9new_phi00[i,] <- fourwise_estimates9[1,]
  simulation_results9new_phi01[i,] <- fourwise_estimates9[3,]
  simulation_results9new_phi10[i,] <- fourwise_estimates9[5,]
  simulation_results9new_phi11[i,] <- fourwise_estimates9[7,]
  
  sim9new_full[,(i+2)] <- ancestor_four_estimation(L = 8, m = 3, four_est = fourwise_estimates9, three_est = threewise_estimates9)
}


## NEIGHBORS
num_nonzero_t <- sum(ancestor_sim8[,1] != 0)
num_nonzero_t

nonzero_bits_t <- which((ancestor_sim8[,1] != 0), arr.ind = TRUE)
zero_bits_t <- which((ancestor_sim8[,1] == 0), arr.ind = TRUE)

neighbors_t <- vector(length = 0) #sequences within a hamming distance of 1 from a true non-zero (95% match)
neighbors2_t <- vector(length = 0) #sequences within a hamming distance of 2 from a true non-zero (90% match)
neighbors3_t <- vector(length = 0) #sequences within a hamming distance of 3 from a true non-zero (85% match)

# FIXED VERSION OF THIS FOR LOOP
# iterate over every possible ancestor sequence
for (i in 1:(dim(an_seq)[1])){
  print(i)
  # calculate hamming distances b/w true non-zeroes and the current sequence
  distances <- t( apply(an_seq[nonzero_bits_t,], 1, FUN =  hamming.distance, y = an_seq[i,]))
  
  # iterate over the distances
  for (j in 1:(length(nonzero_bits_t))){
    # check if distance leq 1
    if (distances[,j] <= 1){
      # add sequence into vector if true
      neighbors_t <- c(neighbors_t,i)
    }
    # check if distance leq 2
    if (distances[,j] <= 2){
      # add sequence into vector if true
      neighbors2_t <- c(neighbors2_t, i) 
    }
    # repeat for leq 3
    if (distances[,j] <= 3){
      # add sequence into vector if true
      neighbors3_t <- c(neighbors3_t, i)
    }
    n3 <- which(distances <= 3, arr.ind = TRUE)
  }
}

# remove repeats from each vector
unique_neighbors_t <- unique(neighbors_t)
unique_neighbors2_t <- unique(neighbors2_t)
unique_neighbors3_t <- unique(neighbors3_t)

# put into increasing order and concatenate with true non-zeroes 
nonzeros_and_neighbors_t <- sort(unique(c(nonzero_bits_t, unique_neighbors_t)))
nonzeros_and_90neighbors_t <- sort(unique(c(nonzero_bits_t, unique_neighbors2_t)))
nonzeros_and_85neighbors_t <- sort(unique(c(nonzero_bits_t, unique_neighbors3_t)))


## CALCULATE DENSITY SUMS
sim9new_full_pf <- matrix(NA, nrow = (2^8), ncol = 6)
colnames(sim9new_full_pf) <- c("Bitwise", "True", "Avg Est", "Bias", "SD", "MSE")

sim9new_full_pf[,(1:2)] <- sim9new_full[,(1:2)]
for (i in 1:(2^8)){
  sim9new_full_pf[i,3] <- mean(sim9new_full[i,(3:102)])
  sim9new_full_pf[i,4] <- sim9new_full_pf[i,3] - sim9_full_pf[i,2] # mean - true
  sim9new_full_pf[i,5] <- sd(sim9new_full[i,(3:102)])
  sim9new_full_pf[i,6] <- sim9new_full_pf[i,5] + sim9new_full_pf[i,4]^2 #sd + bias squared
}

sum(sim9new_full_pf[nonzeros_and_90neighborst, 3])
sum(sim9new_full_pf[nonzeros_and_85neighbors_t, 3])

