df8 <- matrix(NA, nrow = 100, ncol = 5)
colnames(df8) <- c("Rep", "Nonzeros", "HD_1", "HD_2", "HD_3")
df8[,1] <- c(1:100)

for(j in 1:100){
  df8[j,2] <- sum(sim8_full[nonzero_bits, (j+2)])
  df8[j,3] <- sum(sim8_full[nonzeros_and_neighbors, (j+2)])
  df8[j,4] <- sum(sim8_full[nonzeros_and_90neighbors, (j+2)])
  df8[j,5] <- sum(sim8_full[nonzeros_and_85neighbors, (j+2)])
}

plot(x = df8[,1], y = df8[,5], type = "p", xlab = "Replication", ylab = "Total Density Sum", main = "Simulation 8")


df <- matrix(NA, nrow = 100, ncol = 5)
colnames(df) <- c("Rep", "Nonzeros", "HD_1", "HD_2", "HD_3")
df[,1] <- c(1:100)
for (j in 1:100){
  df[j,2] <- sum(sim9_full[nonzero_bits,(j+2)])
  df[j,3] <- sum(sim9_full[nonzeros_and_neighbors,(j+2)])
  df[j,4] <- sum(sim9_full[nonzeros_and_90neighbors,(j+2)])
  df[j,5] <- sum(sim9_full[nonzeros_and_85neighbors, (j+2)])
}

plot(x = df[,1], y = df[,5], type = "p", xlab = "Replication", ylab = "Total Density Sum", main = "Simulation 9")
