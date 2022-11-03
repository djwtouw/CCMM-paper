################################################################################
# Script for analyzing the scalability of CCMM in n.
################################################################################

rm(list = ls())

library(CCMMR)


n = seq(from = 6000, to = 20000, by = 500)
n_reps = 5

# Initialize result
result = matrix(0, nrow = length(n), ncol = 3)
colnames(result) = c("time", "iterations", "n_lambdas")

# Set parameters for the weights
k = 15
phi = 2

path = "Numerical Results/Data/Half Moons/Half_Moons_"

for (i in 1:length(n)) {
    for (j in 1:n_reps) {
        # Read data
        X = as.matrix(read.csv(paste(path, n[i], "/X_", j, ".csv", sep = ""), 
                               header = F))
       
        # Compute weights based on 15 nearest neighbors
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)
        
        # Find lambda for the numbers of clusters around 90% of the data
        lambdas = seq(0, 0.006, 0.000025)
        res = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE,
                                 save_clusterpath = FALSE)
        lambda_idx = which.min(abs(res$info$clusters - 0.9 * n[i]))
        lambda = res$lambdas[lambda_idx]
        
        # Set sequence for lambda
        lambdas = seq(0, lambda, lambda / 100)^2 / lambda
        
        # Compute clusterpath
        res = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE)
        
        # Store the times
        result[i, "time"] = result[i, "time"] + res$elapsed_time
        result[i, "iterations"] = result[i, "iterations"] + sum(res$info[, "iterations"])
        result[i, "n_lambdas"] = result[i, "n_lambdas"] + length(lambdas) - 1
        
        print(result)
    }
}

# Store the results
data = cbind(n, result)
write.table(data, "Numerical Results/Output/scaling_n_ccmm.csv", sep = ",",
            row.names = F)
