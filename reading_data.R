# date created: 06-25-2021
# edited: 06-28-2021
# edited: 07-19-2021

#This script uses population YRI for chromosome 1. See mccl_documentation.Rmd for full documentation

#load packages here:
library(utils)
library(data.table)
library(readr)
library(dplyr)

#######################
# reading in the data:
## unrelated
yri_unr_1 <- read.csv2(file = "C:/Users/gbean/Documents/MCCL REU files/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.gz", header = TRUE, sep = "")
## trios
yri_trio_1 <- read.csv2(file = "C:/Users/gbean/Documents/MCCL REU files/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.phased.gz", header = TRUE, sep = "")
## duos
yri_duo_1 <- read.csv2(file = "C:/Users/gbean/Documents/MCCL REU files/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.D.phased.gz", header = TRUE, sep = "")
#######################

#######################
# SET DATA TO BE USED HERE
hapmap_simulation_data <- yri_trio_1
#######################


######################
# This section codes the full dataset yri_unr_1 into binary sequences, stored in new dataset yri_unr_1_binary

## select everything except for the SNP rsIDs and their locations
snps <- hapmap_simulation_data[,-c(1,2)]
## select the SNP rsIDs and their locations
snps_id_loc <- hapmap_simulation_data[,c(1,2)]

# check that each row has at most 2 letters, not less than 1 letter, and no NA's
unique_letter <- apply(snps, 1, unique)
num_letters <- sapply(unique_letter, length)
print("Num rows w/ invalid number of unique letters:")
sum(num_letters >2 & num_letters < 1)
print("Num rows w/ NA unique letters:")
sum(is.na(num_letters))
print("NOTE both these should be zero!")

## Outline for following code -->
# Step one: count the number of each allele for each SNP
# Step two: identify the major allele
# Step three: check if each allele for the SNP is that allele; if yes, code 0 and code 1 otherwise
# Step four: repeat steps one through three for each row of the dataset and store into new dataset

# Step one/two
## This is a self-written function that takes a SNP row and returns a string of the major allele in that row
### parameter s passed to this function is a SNP, a row of yri_unr_1
find_major_allele <- function(s){
  ## create a vector of the alleles in this row
  alleles_in_snp <- unique(s)
  ## initialize a vector where the length is the number of alleles in the SNP
  num_alleles <- rep(NA, length(alleles_in_snp))
  
  ## calculate the number of each allele in the SNP
  for (j in 1:length(num_alleles)){
    num_alleles[j] <- sum(s %in% alleles_in_snp[j])
  }
  
  ## pick the major allele based on counts
  major_allele <- alleles_in_snp[which.max(num_alleles)]
  return (major_allele)
}

# Step three
## This is a self written function that takes a SNP row and a string of the major allele in that row and returns the row converted
##  to a binary sequence where the major allele is 0 and the minor allele is 1.
### The parameter snp_row is a SNP / row from yri_unr_1
### The parameter allele_m is a string of one letter, ATCG, indicating the major allele for that SNP
binary_seq <- function(snp_row, allele_m){
  seq <- rep(NA, length = length(snp_row))
  for (i in 1:length(snp_row)){
    if (snp_row[i] == allele_m){
      seq[i] <- 0
    }else{
      seq[i] <- 1
    }
  }
  return (seq)
}

# Step four
## Write a for loop iterating over all SNPs in the dataset and calling both functions on each and saving the generated binary sequence 
##  rows to a new dataset

# initialize new dataset as an empty matrix
hapmap_binary <- matrix(NA, nrow = dim(snps)[1], ncol = dim(snps)[2])

# find the major allele for each row
alleles_major <- apply(snps, 1, find_major_allele)

# loop over each SNP
for (i in 1:dim(snps)[1]){
  # generate the binary sequence for that SNP
  binary <- binary_seq(snps[i,], alleles_major[i])
  # attach the recoded row into the initialized dataset
  hapmap_binary[i,] <- binary}

# use the rsIDs from original dataset as the row names
rownames(hapmap_binary) <- snps_id_loc[,1]
# use the genotype IDs from the original dataset as the column names
colnames(hapmap_binary) <- colnames(snps)
######################