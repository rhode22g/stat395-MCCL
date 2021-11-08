# date created: 09-29-2021
# last edited: 11-01-2021
## this script is a new version of a function to estimate the full ancestral distributions using patterns to specify matrices

###############
# LIBRARY CALLS
library(sfsmisc)
library(AlgebraicHaploPackage)


###############
# SET TEST VALUES
# m <- 3
# L <- 15
# test_ancestor <- ancestor(L=L, hapmap_binary)
# test_example <- descendent_sample(L=L, q=0.1, n=5, seed=7, hapmap_binary)
# test_onewise <- estimates_m0(test_example, L=L, n=5)
# test_pairwise <- estimates_m1(L=L, q=0.1, n=5, d=test_example)
# test_threewise <- estimates_m2(L=L, q=0.1, n=5, descend = test_example)
# test_fourwise <- estimates_m3(L = L, q = 0.1, n=5, d = test_example)
###############


#############
# SCRATCH WORK FOR PAIRWISE
# creating matrices

# pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
# for (i in 1:(L-m)){
#   pairs_needed[,i] <- c(rep(test_pairwise[1,i], times = (1/2^(m+i))*2^L), 
#                         rep(test_pairwise[2,i], times = (1/2^(m+i))*2^L), 
#                         rep(test_pairwise[3,i], times = (1/2^(m+i))*2^L), 
#                         rep(test_pairwise[4,i], times = (1/2^(m+i))*2^L))
# }
# 
# 
# 
# ones_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
# for (j in 1:(L-m-1)){
#   ones_needed[,j] <- rep(c(rep(test_onewise[1,(j+1)], times = (1/2^(m+j))*2^L),
#                            rep(test_onewise[2,(j+1)], times = (1/2^(m+j))*2^L)), times = 2)
# }
# 
# # ones_needed <- as.matrix(ones_needed)
# 
# products <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
# products[,1:(L-m)] <- pairs_needed
# products[,(L-m+1):(2*L - 2*m - 1)] <- (1/ones_needed)
# 
# ancestor_estimates <- matrix(NA, nrow = (2^L), ncol = 1)
# ancestor_estimates <- apply(products, 1, prod)
# ancestor_estimates[is.na(ancestor_estimates)] <- 0
# ancestor_estimates <- as.matrix(ancestor_estimates)
# 
# sum(ancestor_estimates)


##########
# FROM PAIRWISE
#########
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


##########
# FROM THREEWISE
#########
ancestor_three_estimation <- function(L, m, three_est, pairs_est){
  
  # define matrix of needed threewise estimates
  threes_needed <- matrix(NA, nrow = (2^L) , ncol = (L-m))
  for (i in 1:(L-m)){
    threes_needed[,i] <- c(rep(three_est[1,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[2,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[3,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[4,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[5,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[6,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[7,i], times = (1/2^(m+i))*2^L),
                           rep(three_est[8,i], times = (1/2^(m+i))*2^L))
  }
  
  # define matrix of needed pairwise estimates
  pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  for(j in 1:(L-m-1)){
    pairs_needed[,j] <- c(rep(pairs_est[1,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[2,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[3,(j+1)], times = (1/2^(m+j))*2^L),
                          rep(pairs_est[4,(j+1)], times = (1/2^(m+j))*2^L))
  }
  
  # concatenate into one matrix and invert pairs
  to_mult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m -1))
  to_mult[,(1:(L-m))] <- threes_needed
  to_mult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(pairs_needed)
  
  # calculate estimates
  an_est <- matrix(NA, nrow = (2^L), ncol = 1)
  an_est <- apply(to_mult, 1, prod)
  
  # replace NA's with 0
  an_est[is.na(an_est)] <- 0
  
  # return matrix of elements
  an_est <- as.matrix(an_est)
  return(an_est)
}

#############
# SCRATCH WORK THREEWISE
# threes_needed <- matrix(NA, nrow = (2^L) , ncol = (L-m))
# for (i in 1:(L-m)){
#   threes_needed[,i] <- c(rep(test_threewise[1,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[2,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[3,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[4,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[5,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[6,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[7,i], times = (1/2^(m+i))*2^L),
#                          rep(test_threewise[8,i], times = (1/2^(m+i))*2^L))
# }
# 
# pairs_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
# for(j in 1:(L-m-1)){
#   pairs_needed[,j] <- c(rep(test_pairwise[1,(j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_pairwise[2,(j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_pairwise[3,(j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_pairwise[4,(j+1)], times = (1/2^(m+j))*2^L))
# }
# 
# to_mult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m -1))
# to_mult[,(1:(L-m))] <- threes_needed
# to_mult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(pairs_needed)
# 
# an_est <- matrix(NA, nrow = (2^L), ncol = 1)
# an_est <- apply(to_mult, 1, prod)
# an_est[is.na(an_est)] <- 0
# an_est <- as.matrix(an_est)

###################


#####################
# FROM FOURWISE ####
####################
ancestor_four_estimation <- function(L, m, four_est, three_est){
  # initialize matrix for needed fourwise estimates
  fours_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
  
  # fill matrix with needed fourwise estimates
  for (i in 1:(L-m)){
    fours_needed[,i] <- c(rep(four_est[1,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[2,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[3,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[4,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[5,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[6,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[7,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[8,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[9,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[10,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[11,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[12,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[13,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[14,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[15,i], times = (1/2^(m+i))*2^L),
                         rep(four_est[16,i], times = (1/2^(m+i))*2^L))
  }
  
  # initialize matrix for needed threewise estimates
  threes_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
  
  # fill matrix with needed threewise estimates
  for (j in 1:(L-m-1)){
    threes_needed[,j] <- c(rep(three_est[1, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[2, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[3, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[4, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[5, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[6, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[7, (j+1)], times = (1/2^(m+j))*2^L),
                          rep(three_est[8, (j+1)], times = (1/2^(m+j))*2^L))
  }
  
  # concatenate into one matrix & take inverse of threewise
  to_mult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
  to_mult[,(1:(L-m))] <- fours_needed
  to_mult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(threes_needed)
  
  # calculate estimates
  an_ests <- matrix(NA, nrow = (2^L), ncol = 1)
  an_ests <- apply(to_mult, 1, prod)
  
  # replace NA's with zeros
  an_ests[is.na(an_ests)] <- 0
  
  # return matrix of estimates
  an_ests <- as.matrix(an_ests)
  return(an_ests)
}
#####################


####################
# SCRATCH WORK FOR FOURWISE
# four_needed <- matrix(NA, nrow = (2^L), ncol = (L-m))
# for (i in 1:(L-m)){
#   four_needed[,i] <- c(rep(test_fourwise[1,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[2,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[3,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[4,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[5,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[6,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[7,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[8,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[9,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[10,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[11,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[12,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[13,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[14,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[15,i], times = (1/2^(m+i))*2^L),
#                        rep(test_fourwise[16,i], times = (1/2^(m+i))*2^L))
# }
# 
# three_needed <- matrix(NA, nrow = (2^L), ncol = (L-m-1))
# for (j in 1:(L-m-1)){
#   three_needed[,j] <- c(rep(test_threewise[1, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[2, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[3, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[4, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[5, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[6, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[7, (j+1)], times = (1/2^(m+j))*2^L),
#                         rep(test_threewise[8, (j+1)], times = (1/2^(m+j))*2^L))
# }
# 
# tomult <- matrix(NA, nrow = (2^L), ncol = (2*L - 2*m - 1))
# tomult[,(1:(L-m))] <- four_needed
# tomult[,((L-m+1):(2*L - 2*m - 1))] <- 1/(three_needed)
# 
# estimates <- matrix(NA, nrow = (2^L), ncol = 1)
# estimates <- apply(tomult, 1, prod)
# estimates[is.na(estimates)] <- 0
# 
# estimates <- as.matrix(estimates)
#######################

####################
# test function performance
# ancestor_test <- ancestor_pair_estimation(L=L, m=m, pairs_est = test_pairwise, ones_est = test_onewise)
# ancestor_test_2 <- ancestor_three_estimation(L = L, m=2, three_est = test_threewise, pairs_est = test_pairwise)
# ancestor_test_3 <- ancestor_four_estimation(L = L, m = 3, four_est = test_fourwise, three_est = test_threewise)
###################