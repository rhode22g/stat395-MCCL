# date created: 09-29-2021
# last edited: 10-04-2021
## this script is a new version of a function to estimate the full ancestral distributions using patterns to specify matrices

#test values
library(sfsmisc)
library(AlgebraicHaploPackage)

m <- 1
L <- 20
test_ancestor <- ancestor(L=L, hapmap_binary)
test_example <- descendent_sample(L=L, q=0.1, n=5, seed=7, hapmap_binary)
test_onewise <- estimates_m0(test_example, L=L, n=5)
test_pairwise <- estimates_m1(L=L, q=0.1, n=5, d=test_example)


# creating matrices
# pairs_needed <- data.frame(site_1_2 = c(rep(test_pairwise[1,1], times = (1/2^(m+1))*2^L), 
#                              rep(test_pairwise[2,1], times = (1/2^(m+1))*2^L), 
#                              rep(test_pairwise[3,1], times = (1/2^(m+1))*2^L), 
#                              rep(test_pairwise[4,1], times = (1/2^(m+1))*2^L)),
#                            
#                            site_2_3 = c(rep(test_pairwise[1,2], times = (1/2^(m+2))*2^L), 
#                              rep(test_pairwise[2,2], times = (1/2^(m+2))*2^L), 
#                              rep(test_pairwise[3,2], times = (1/2^(m+2))*2^L), 
#                              rep(test_pairwise[4,2], times = (1/2^(m+2))*2^L)),
#                            
#                            site_3_4 = c(rep(test_pairwise[1,3], times = (1/2^(m+3))*2^L), 
#                              rep(test_pairwise[2,3], times = (1/2^(m+3))*2^L), 
#                              rep(test_pairwise[3,3], times = (1/2^(m+3))*2^L), 
#                              rep(test_pairwise[4,3], times = (1/2^(m+3))*2^L)),
#                            
#                            site_4_5 = c(rep(test_pairwise[1,4], times = (1/2^(m+4))*2^L), 
#                              rep(test_pairwise[2,4], times = (1/2^(m+4))*2^L), 
#                              rep(test_pairwise[3,4], times = (1/2^(m+4))*2^L), 
#                              rep(test_pairwise[4,4], times = (1/2^(m+4))*2^L)))

pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
for (i in 1:(L-m)){
  pairs_needed[,i] <- c(rep(test_pairwise[1,i], times = (1/2^(m+i))*2^L), 
                        rep(test_pairwise[2,i], times = (1/2^(m+i))*2^L), 
                        rep(test_pairwise[3,i], times = (1/2^(m+i))*2^L), 
                        rep(test_pairwise[4,i], times = (1/2^(m+i))*2^L))
}

# pairs_needed <- as.matrix(pairs_needed)

# ones_needed <- data.frame(site_1 = rep(c(rep(test_onewise[1,2], times = (1/2^(m+1))*2^L),
#                                      rep(test_onewise[2,2], times = (1/2^(m+1))*2^L)), times = 2),
#                           
#                           site_2 = c(rep(test_onewise[1,3], times = (1/2^(m+2))*2^L),
#                                      rep(test_onewise[2,3], times = (1/2^(m+2))*2^L)),
#                           
#                           site_3 = c(rep(test_onewise[1,4], times = (1/2^(m+3))*2^L),
#                                      rep(test_onewise[2,4], times = (1/2^(m+3))*2^L)))

ones_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
for (j in 1:(L-m-1)){
  ones_needed[,j] <- rep(c(rep(test_onewise[1,(j+1)], times = (1/2^(m+j))*2^L),
                           rep(test_onewise[2,(j+1)], times = (1/2^(m+j))*2^L)), times = 2)
}

# ones_needed <- as.matrix(ones_needed)

products <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
products[,1:(L-m)] <- pairs_needed
products[,(L-m+1):(2*L - 2*m - 1)] <- (1/ones_needed)

ancestor_estimates <- matrix(NA, nrow = (2^L), ncol = 1)
ancestor_estimates <- apply(products, 1, prod)
ancestor_estimates[is.na(ancestor_estimates)] <- 0
ancestor_estimates <- as.matrix(ancestor_estimates)

sum(ancestor_estimates)


ancestor_pair_estimation <- function(L, m, pairs_est, ones_est){
  # define matrix of needed pair estimates
  pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
  for (i in 1:(L-m)){
    pairs_needed[,i] <- c(rep(pairs_est[1,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[2,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[3,i], times = (1/2^(m+i))*2^L), 
                          rep(pairs_est[4,i], times = (1/2^(m+i))*2^L))
  }
  
  # define matrix of needed onewise estimates
  ones_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  for (j in 1:(L-m-1)){
    ones_needed[,j] <- rep(c(rep(ones_est[1,(j+1)], times = (1/2^(m+j))*2^L),
                             rep(ones_est[2,(j+1)], times = (1/2^(m+j))*2^L)), times = 2)
  }
  
  # concatenate into one matrix
  products <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
  products[,1:(L-m)] <- pairs_needed
  products[,(L-m+1):(2*L - 2*m - 1)] <- (1/ones_needed)
  
  # calculate estimates and replace NA's with 0
  ancestor_estimates <- matrix(NA, nrow = (2^L), ncol = 1)
  ancestor_estimates <- apply(products, 1, prod)
  ancestor_estimates[is.na(ancestor_estimates)] <- 0
  ancestor_estimates <- as.matrix(ancestor_estimates)
  
  # return matrix of estimates
  return(ancestor_estimates)
}

# test function performance
ancestor_test <- ancestor_pair_estimation(L=L, m=m, pairs_est = test_pairwise, ones_est = test_onewise)
