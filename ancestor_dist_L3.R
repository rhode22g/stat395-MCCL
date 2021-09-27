# date created: 06-29-2021
## This script obtains the ancestral distribution for L=3 using the HapMap data as the ancestor pool. See mccl_documentation.Rmd for full
##  documentation.

###############
# loading packages
library(dplyr)
###############

###############
# We know that the possible sequences for L=3 are the following 8 sequences
#  (0,0,0) , (0,0,1) , (0,1,1) , (0,1,0), (1,0,0), (1,0,1), (1,1,1), (1,1,0)
# Assign these sequence numbers -->
### 1: (0,0,0)
### 2: (0,0,1)
### 3: (0,1,0)
### 4: (0,1,1)
### 5: (1,0,0)
### 6: (1,0,1)
### 7: (1,1,0)
### 8: (1,1,1)
###############

###############
# Take the subset of the binary unrelated dataset so that each haplotype is L=3
unr_bin_L3 <- yri_unr_1_binary[1:3,]
###############

###############
# Create a new row that stores each length 3 haplotype as a single string

# initialize empty row
hap_3 <- rep(NA, length = dim(unr_bin_L3)[1])
# combine new row with dataset of first three SNPs
new_unr_bin_L3 <- rbind(unr_bin_L3, hap_3)

# iterate over COLUMNS
for (i in 1:dim(unr_bin_L3)[2]){
  hap <- vector(length = 3)
  # iterate over ROWS
  for (j in 1:dim(unr_bin_L3)[1]){
    hap[j] <- unr_bin_L3[j,i]
    # print(hap)
    # print(unr_bin_L3[j,i])
  }
  full_hap <- toString(hap)
  new_unr_bin_L3[4,i] <- full_hap}
###############

###############
# Calculate proportions

## set sequences as strings and initialize counts of each sequence at zero
# s1 <- "0, 0, 0"
# num1 <- 0
# s2 <- "0, 0, 1"
# num2 <- 0
# s3 <- "0, 1, 0"
# num3 <- 0
# s4 <- "0, 1, 1"
# num4 <- 0
# s5 <- "1, 0, 0"
# num5 <- 0
# s6 <- "1, 0, 1"
# num6 <- 0
# s7 <- "1, 1, 0"
# num7 <- 0
# s8 <- "1, 1, 1"
# num8 <- 0

## instead of defining each variable by hand, use iteration to create a vector of strings
sequences <- vector(length = 8)
i <- 1
for (k in 0:1){
  for (h in 0:1){
    for (q in 0:1){
      sequences[i] <- toString(c(k,h,q))
      i <- i + 1
    }
  }
}

## instead of initializing each count by hand, use iteration to create a vector of counts initialized as 0's
sequence_counts <- vector(length = 8)
for (v in 1:8){
  sequence_counts[v] <- 0
}

## pull out the row of each haplotype of the format "x1, x2, x3"
haps_length3 <- new_unr_bin_L3[4,]

## iterate over the haplotype strings, check if they match the sequence variables and increment counts accordingly
# for (x in 1:length(haps_length3)){
#   if (haps_length3[x] == s1){
#     num1 <- num1 + 1
#   }
#   if (haps_length3[x] == s2){
#     num2 <- num2 + 1
#   }
#   if (haps_length3[x] == s3){
#     num3 <- num3 + 1
#   }
#   if (haps_length3[x] == s4){
#     num4 <- num4 + 1
#   }
#   if (haps_length3[x] == s5){
#     num5 <- num5 + 1
#   }
#   if (haps_length3[x] == s6){
#     num6 <- num6 + 1
#   }
#   if (haps_length3[x] == s7){
#     num7 <- num7 + 1
#   }
#   if (haps_length3[x] == s8){
#     num8 <- num8 + 1
#   }
# }

## alternative to separate if statements
### iterate over haplotype strings
for (x in 1:length(haps_length3)){
  for (s in 1:length(sequences)){
    if (haps_length3[x] == sequences[s]){
      sequence_counts[s] <- sequence_counts[s] + 1
    }
  }
}

## create a a vector of the counts for each sequence
# sequence_counts <- c(num1, num2, num3, num4, num5, num6, num7, num8)


## initialize an empty vector to hold the sequence proportions
sequence_proportions <- vector(length = 8)

## iterate over each sequence count and calculate the proportion for each sequence as the seq count divided by num of haplotypes
for (s in 1:length(sequence_counts)){
  sequence_proportions[s] <- sequence_counts[s]/length(haps_length3)}

###############