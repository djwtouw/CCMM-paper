################################################################################
# Script that tracks the computation for AMA and CCMM. Clusterpaths are desgined
# such that no fusions occur and the comparison of computation times focuses on
# the efficiency of the update.
################################################################################

rm(list = ls())

# Load libraries, be sure to restart R each time when cvxclustr is reloaded
library(CCMMR)
library(cvxclustr)

# Extra functions for converting AMA result
source("src/R/convert_to_cp.R")
source("src/R/save_load.R")
Rcpp::sourceCpp("src/Cpp/utils.cpp")

# Data sizes
N = seq(500, 20000, 500)
n_reps = 10
n_lambdas = 5

# Set parameters for the weights
k = 15
phi = 0.5

# Generate and store data
set.seed(1)
for (n_i in 1:length(N)) {
    for (r in 1:n_reps) {
        # Number of observations
        n = N[n_i]
        
        # File name
        fname = paste("Numerical Results/Data/Grid/X_", n, "_", r, ".csv", 
                      sep = "")
        
        # If the file exists already, skip to the next iteration
        if (file.exists(fname)) {
            next
        }
        
        # Dimension 1 of the grid
        p = 20
        
        # Create the data
        X = matrix(0, nrow = n, ncol = 2)
        for (i in 1:(n / p)) {
            for (j in 1:p) {
                X[(i - 1) * p + j, 1] = i + rnorm(1) / 15
                X[(i - 1) * p + j, 2] = j + rnorm(1) / 15
            }
        }
        
        # Write the data to disc
        write.table(X, file = fname, sep = ",", col.names = F, row.names = F)
    }
}

# For each data set, find a value for lambda for which the number of clusters is
# reduced by one. Store a value that is just below that value that can be used
# in the comparisons of the algorithms
for (n_i in 1:length(N)) {
    for (r in 1:n_reps) {
        # Load X
        fname = paste("Numerical Results/Data/Grid/X_", N[n_i], "_", r, ".csv", 
                      sep = "")
        X = as.matrix(read.csv(fname, header = FALSE))
        
        # Get sparse distances in dictionary of keys format
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)

        # Find the smallest lambda for which n - 1 clusters remain
        ccmm = convex_clustering(X, W, N[n_i] - 5, N[n_i], lambda_init = 0.001, 
                                 scale = FALSE, center = FALSE,
                                 max_iter_phase_2 = 50)
        
        # The lambda chosen for this instance should be smaller than the value
        # for which fewer than n clusters remained
        lambda = 0.9 * ccmm$lambdas[2]
        
        # Print lambda
        cat(paste("Lambda for n = ", N[n_i], "(" , r, ") is ", 
                  round(lambda, 3), "\n", sep = ""))
        
        # File name
        fname = paste("Numerical Results/Data/Grid/lambda_", N[n_i], "_", r,
                      ".bin", sep = "")
        
        # If the file exists already, check if the value for lambda is the same
        if (file.exists(fname)) {
            if (abs(as.numeric(load_bin(fname)) - lambda) < 1e-10) {
                next
            }
        }
        
        # Write this value for lambda to disc
        save_bin(seq(0, lambda, length.out = n_lambdas), fname)
    }
}

################################################################################
# Now run the code in no_clustering_timings.m
################################################################################

################################################################################
# CCMM vs SSNAL
################################################################################
# Create matrix to store computation times
timings_ccmm_ssnal = matrix(0, nrow = length(N), ncol = 2)
colnames(timings_ccmm_ssnal) = c("CCMM", "SSNAL")

# Create matrices to store the loss function values
losses_ccmm_ssnal = matrix(0, nrow = length(N), ncol = 2)
colnames(losses_ccmm_ssnal) = c("CCMM", "SSNAL")

# SSNAL Timings
timings_ssnal = read.csv(paste("Numerical Results/Output/no_clustering_timings",
                               "_ssnal.csv", sep = ""), header = FALSE)
timings_ccmm_ssnal[, "SSNAL"] = timings_ssnal[, 1]
rm(timings_ssnal)

# CCMM vs SSNAL
for (n_i in 1:length(N)) {
    # Iterate over the different realizations of the data for each size
    for (r in 1:n_reps) {
        # Read in the data
        fname = paste("Numerical Results/Data/Grid/X_", N[n_i], "_", r, ".csv", 
                      sep = "")
        X = as.matrix(read.csv(fname, header = FALSE))
        
        # Read in the value for lambda
        fname = paste("Numerical Results/Data/Grid/lambda_", N[n_i], "_", r,
                      ".bin", sep = "")
        lambdas = c(load_bin(fname, n = n_lambdas))
        
        # Read in losses attained by SSNAL
        fname = paste("Numerical Results/Output/Grid/losses_", N[n_i], "_", r,
                      ".bin", sep = "")
        losses_ssnal = c(load_bin(fname, n = n_lambdas))
        
        # Get sparse distances in dictionary of keys format
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)
        
        # Compute results CCMM
        ccmm = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE,
                                  target_losses = losses_ssnal)
        timings_ccmm_ssnal[n_i, "CCMM"] = timings_ccmm_ssnal[n_i, "CCMM"] +
            ccmm$elapsed_time
        
        # Store the loss values in the matrix
        losses_ccmm_ssnal[n_i, "CCMM"] = losses_ccmm_ssnal[n_i, "CCMM"] +
            ccmm$info[n_lambdas, "loss"]
        losses_ccmm_ssnal[n_i, "SSNAL"] = losses_ccmm_ssnal[n_i, "SSNAL"] +
            losses_ssnal[n_lambdas]
        
        # Print timings to keep track of progress
        print(timings_ccmm_ssnal)
    }
}

# Calculate the ratio of SSNAL to CCMM
ratio = losses_ccmm_ssnal[, "SSNAL"] / losses_ccmm_ssnal[, "CCMM"]

# Create the scatter plot
og = par(mfrow = c(1, 2))
plot(N, ratio, pch = 16, cex = 1, col = "blue", xlab = "N",
     ylab = "SSNAL / CCMM Ratio", main = "SSNAL / CCMM Ratio (Losses)", 
     ylim = c(min(min(ratio), 1) * (1 - 1e-6), max(max(ratio), 1) * (1 + 1e-6)))

# Add a reference line at ratio = 1
abline(h = 1, col = "red", lty = 2)

# Calculate the ratio of SSNAL to CCMM
ratio = timings_ccmm_ssnal[, "SSNAL"] / timings_ccmm_ssnal[, "CCMM"]

# Create the scatter plot
plot(N, ratio, pch = 16, cex = 1, col = "blue", xlab = "N",
     ylab = "SSNAL / CCMM Ratio", main = "SSNAL / CCMM Ratio (Timings)", 
     ylim = c(0, max(max(ratio), 1) * 1.1))
par(og)

write.table(
    timings_ccmm_ssnal,
    "Numerical Results/Output/no_clustering_timings_ccmm_ssnal.csv",
    row.names = FALSE, col.names = TRUE, sep = ","
)

write.table(
    losses_ccmm_ssnal,
    "Numerical Results/Output/no_clustering_losses_ccmm_ssnal.csv",
    row.names = FALSE, col.names = TRUE, sep = ","
)

################################################################################
# CCMM versus AMA
################################################################################
# Create matrix to store computation times
timings_ccmm_ama = matrix(0, nrow = length(N), ncol = 2)
colnames(timings_ccmm_ama) = c("CCMM", "AMA")

# Create matrices to store the loss function values
losses_ccmm_ama = matrix(0, nrow = length(N), ncol = 2)
colnames(losses_ccmm_ama) = c("CCMM", "AMA")

# CCMM vs AMA
for (n_i in 1:19) {
    # Iterate over the different realizations of the data for each size
    for (r in 1:n_reps) {
        # Read in the data
        fname = paste("Numerical Results/Data/Grid/X_", N[n_i], "_", r, ".csv", 
                      sep = "")
        X = as.matrix(read.csv(fname, header = FALSE))
        
        # Read in the value for lambda
        fname = paste("Numerical Results/Data/Grid/lambda_", N[n_i], "_", r,
                      ".bin", sep = "")
        lambdas = c(load_bin(fname, n = n_lambdas))
        
        ########################################################################
        # AMA
        ########################################################################
        # Get the weights in vector format
        w_vec = kernel_weights(t(X), phi)
        w_vec = knn_weights(w_vec, k, nrow(X))
        
        # Compute recommended step size nu
        nu = AMA_step_size(w_vec, nrow(X))
        
        # Compute AMA result
        t1 = Sys.time()
        ama = cvxclust_path_ama(t(X), w_vec, lambdas, nu = nu)
        t2 = Sys.time()
        timings_ccmm_ama[n_i, "AMA"] = timings_ccmm_ama[n_i, "AMA"] +
            as.numeric(difftime(t2, t1, units = "secs"))
        
        # Compute clusterpath from AMA result
        ama_cp = convert_to_cp(ama, lambdas)
        losses_ama = loss_from_cp(ama_cp$coordinates, X, w_vec, lambdas)[, 1]
        
        # Store the loss values in the matrix
        losses_ccmm_ama[n_i, "AMA"] = losses_ccmm_ama[n_i, "AMA"] +
            losses_ama[n_lambdas]
        
        ########################################################################
        # CCMM
        ########################################################################
        # Get sparse distances in dictionary of keys format
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)
        
        # Compute results CCMM
        ccmm = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE,
                                  target_losses = losses_ama)
        timings_ccmm_ama[n_i, "CCMM"] = timings_ccmm_ama[n_i, "CCMM"] +
            ccmm$elapsed_time
        
        # Store the loss values in the matrix
        losses_ccmm_ama[n_i, "CCMM"] = losses_ccmm_ama[n_i, "CCMM"] +
            ccmm$info[n_lambdas, "loss"]
        
        # Print timings to keep track of progress
        print(timings_ccmm_ama)
    }
    
    # Make sure to save each iteration
    write.table(
        timings_ccmm_ama,
        "Numerical Results/Output/no_clustering_timings_ccmm_ama.csv",
        row.names = FALSE, col.names = TRUE, sep = ","
    )
    
    write.table(
        losses_ccmm_ama,
        "Numerical Results/Output/no_clustering_losses_ccmm_ama.csv",
        row.names = FALSE, col.names = TRUE, sep = ","
    )
}

# Load results
timings_ccmm_ama = as.matrix(
    read.csv("Numerical Results/Output/no_clustering_timings_ccmm_ama.csv")
)
losses_ccmm_ama = as.matrix(
    read.csv("Numerical Results/Output/no_clustering_losses_ccmm_ama.csv")
)

# Adjust for obtained results
N = N[timings_ccmm_ama[, "AMA"] > 0]
losses_ccmm_ama = losses_ccmm_ama[timings_ccmm_ama[, "AMA"] > 0, ]
timings_ccmm_ama = timings_ccmm_ama[timings_ccmm_ama[, "AMA"] > 0, ]

# Calculate the ratios for the losses
ratio = losses_ccmm_ama[, "AMA"] / losses_ccmm_ama[, "CCMM"]

# Create the scatter plot
og = par(mfrow = c(1, 2))
plot(N, ratio, pch = 16, cex = 1, col = "blue", xlab = "N",
     ylab = "AMA / CCMM Ratio", main = "AMA / CCMM Ratio (Losses)", 
     ylim = c(min(min(ratio), 1) * (1 - 1e-6), max(max(ratio), 1) * (1 + 1e-6)))

# Add a reference line at ratio = 1
abline(h = 1, col = "red", lty = 2)

# Calculate the ratio of AMA to CCMM
ratio = timings_ccmm_ama[, "AMA"] / timings_ccmm_ama[, "CCMM"]

# Create the scatter plot
plot(N, ratio, pch = 16, cex = 1, col = "blue", xlab = "N",
     ylab = "AMA / CCMM Ratio", main = "AMA / CCMM Ratio (Timings)", 
     ylim = c(0, max(max(ratio), 1) * 1.1))
par(og)
