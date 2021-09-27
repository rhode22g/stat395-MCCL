# date created: Aug 15, 2021

# this script defines functions to calculate full ancestral distribution estimates from the (m+1)-wise marginal estimates
library(sfsmisc)
########
# FROM PAIRWISE

# function inputs: 
## length of the sequence
## matrix of pairwise estimates
## matrix of onewise estimates

ancestor_pairs <- function(L, pairs, ones){
  # create a matrix of the ancestor sequences we need estimates of
  an_seq <- matrix(c(0:(2^(L)-1)), ncol=1)
  an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=L) ) 
  an_seq <- matrix( as.character(an_seq), nrow=2^{L}, ncol=L)
  
  # initialize matrix to hold final estimates
  estimates <- matrix(NA, nrow=dim(an_seq)[1], ncol = 1)
  
  # iterate over the rows
  for (i in 1:(dim(an_seq)[1])){
    # initialize vector to hold pair estimates to multiply
    pi_pair_est <- vector(length=(L-1))
    # initialize vector to hold ones estimates to multiply
    pi_ones_est <- vector(length=(L-2))
    # iterate over the columns
    # take needed pair estimates
    for (j in 1:(L-1)){
      # calculate bitwise value for the pair
      site1 <- bitwise(as.numeric(an_seq[i,j]),(2-1))
      site2 <- bitwise(as.numeric(an_seq[i,(j+1)]),(2-2))
      pairs_val <- site1 + site2
      # print(pairs_val)
      pi_pair_est[j] <- pairs[(pairs_val + 1),j]
    }
    # take needed onewise estimates
    for (k in 2:(L-1)){
      site <- as.numeric(an_seq[i,k])
      pi_ones_est[(k-1)] <- ones[(site + 1),k]
    }
    # calculate the estimate
    if (prod(pi_ones_est) == 0){
      estimates[i,] <- 0
    }
    else {
      estimates[i,] <- prod(pi_pair_est) / prod(pi_ones_est)  
    }
  }
  return(estimates)
}
########

#####
# version of above function that calculates only the true non-zeroes
## L is the sequence length, pairs is pairwise estimates, ones is the onewise estimates
## true_an is the matrix of the simulated ancestor probabilities, seq_an is a matrix of the possible sequences for L
ancestor_pairs_nonzero <- function(L, pairs, ones, true_an, seq_an){
  # find the number of true nonzeroes
  numnonzero <- 2^(L) - sum(true_an[,1] == 0)
  
  # initialize matrix to store the estimates
  estimates <- matrix(NA, nrow = numnonzero, ncol = 1)
  
  for (i in 1:dim(seq_an)[1]){
    if (seq_an[i,] != 0){
      
    }
  }
}

#######
# TEST VALUES
test_an_pairs <- ancestor_pairs(L=4, pairs = test_pairwise, ones = test_onewise)

sim_an_results <- matrix(NA, nrow = 100, ncol = 2^4)
sim_an_results[1,] <- test_an_pairs[,1]
######

# an_seq <- matrix(c(0:(2^(4)-1)), ncol=1)
# an_seq <- t( apply(an_seq,1,digitsBase, base=2, ndigits=4) )
# an_seq <- matrix( as.character(an_seq), nrow=2^{4}, ncol=4)
# # an_seq <- apply(an_seq, 1, paste, collapse="")
# 
# L <- 4
# estimates <- matrix(NA, nrow=dim(an_seq)[1], ncol = 1)
# # iterate over the rows
# for (i in 1:(dim(an_seq)[1])){
#   # print("start row:")
#   pi_pair_est <- vector(length=(L-1))
#   pi_ones_est <- vector(length=(L-2))
#   # iterate over the columns
#   # take needed pair estimates
#   for (j in 1:(L-1)){
#     # calculate bitwise value for the pair
#     site1 <- bitwise(as.numeric(an_seq[i,j]),(2-1))
#     site2 <- bitwise(as.numeric(an_seq[i,(j+1)]),(2-2))
#     pairs_val <- site1 + site2
#     # print(pairs_val)
#     pi_pair_est[j] <- test_pairwise[(pairs_val + 1),j]
#   }
#   # take needed onewise estimates
#   for (k in 2:(L-1)){
#     site <- as.numeric(an_seq[i,k])
#     print(site)
#     pi_ones_est[(k-1)] <- test_onewise[(site + 1),k]
#   }
#   # calculate the estimate
#   print(pi_pair_est)
#   print(pi_ones_est)
#   estimates[i,] <- prod(pi_pair_est) / prod(pi_ones_est)
# }
########
# FROM THREEWISE
########
