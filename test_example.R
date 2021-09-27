# date created: July 7, 2021
# date edited: July 16, 2021

# this script finds the m=0,1,2,3 (m+1)-wise marginal estimates for descendant matrix, "test_example," where the ancestor distribution is the YRI
# unrelated population on chromosome 1, q=0.1, L=4, and n=5. the seed is set at 7.

#######
library(sfsmisc)
library(AlgebraicHaploPackage)
test_ancestor <- ancestor(L=5, hapmap_binary)
test_example <- descendent_sample(L=5, q=0.1, n=5, seed=7, hapmap_binary)
test_onewise <- estimates_m0(test_example, L=5, n=5)
test_pairwise <- estimates_m1(L=5, q=0.1, n=5, d=test_example)
test_threewise <- estimates_m2(L=4, q=0.1, n=5, descend = test_example)
test_fourwise <- estimates_m3(L=4, q=0.1, n=5, d=test_example)

test_an_seq <- matrix(c(0:(2^(5)-1)), ncol=1)
test_an_seq <- t( apply(test_an_seq,1,digitsBase, base=2, ndigits=5) )
#######

# A1_0 <- A_1(c(2,3,0, 0), c(0.81, -0.81, -0.81, 0.81))
# A2_0 <- A_2(c(2,3,0, 0), c(0.81, -0.81, -0.81, 0.81), c(0.076, 0.924, 0.324, -0.324))
# A3_0 <- A_3(c(2,3,0, 0), c(0.81, -0.81, -0.81, 0.81), c(0.076, 0.924, 0.324, -0.324))
# A4_0 <- A_4(c(2,3,0, 0), c(0.81, -0.81, -0.81, 0.81), c(0.076, 0.924, 0.324, -0.324))
A1_0 <- 2.152336
A2_0 <- -2.582803
A3_0 <- 1.033121
A4_0 <- -0.1377495
print(A1_0)
print(A2_0)
print(A3_0)
print(A4_0)

# calculate the roots
all_roots <- polyroot(c(A4_0, A3_0, A2_0, A1_0))
print("phi0 all roots")
print(all_roots)

# find the real roots 
real_roots <- pick_real(all_roots)
print("phi0 real roots")
print(real_roots)

test_matrix <- rep(NA, 2^100)

## testing NA and negatives that showed up in simulation
descend_test <- descendent_sample(L=20, q=0.1, n=100, seed =35, hapmap = hapmap_binary)
onewise_test <- estimates_m0(descendents = descend_test, L=20, n=100)
pairwise_test <- estimates_m1(L=20, q=0.1, n=100, d= descend_test)
threewise_test <- estimates_m2(L=20, q=0.1, n=100, descend = descend_test)


# TEST RESOLVED
descend_test_81 <- descendent_sample(L=20, q=0.005, n=100, seed = 81, hapmap_binary)
threewisetest_81 <- estimates_m2(L=20, q=0.005, n=100, descend = descend_test_81)
fourwisetest_81 <- estimates_m3(L=20, q=0.005, n=100, d = descend_test_81)

