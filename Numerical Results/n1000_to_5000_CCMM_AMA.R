################################################################################
# Script that tracks the computation time for CCMM and AMA for n = 1000 to 
# n = 5000. For each n there are 10 data sets. At the end of the script the 
# results obtained for SSNAL are added.
#
# Note: run this script after running [Numerical Results/n1000_to_5000_SSNAL.m] 
#       to ensure that SSNAL results are also processed.
################################################################################

rm(list = ls())

# Load libraries, be sure to restart R each time when cvxclustr is reloaded
library(CCMMR)
library(cvxclustr)

# Extra functions for converting AMA result
source("src/R/convert_to_cp.R")
Rcpp::sourceCpp("src/Cpp/utils.cpp")


# Set parameters for the weights
k = 15
phi = 2
lambdas = seq(0, 110, 0.2)

# Select the sizes of the data
N = seq(1000, 5000, 1000)

# Create matrix to store computation times
timings = matrix(0, nrow = length(N), ncol = 4)
colnames(timings) = c("CCMM", "AMA", "SSNAL", "completed")

# Create matrices to store the loss function values
losses1000 = matrix(0, nrow = length(lambdas), ncol = 3)
losses5000 = matrix(0, nrow = length(lambdas), ncol = 3)
colnames(losses1000) = c("CCMM", "AMA", "SSNAL")
colnames(losses5000) = c("CCMM", "AMA", "SSNAL")

# The base path to the folder where the data is stored
data_path_base = "Numerical Results/Data/Half Moons/Half_Moons_"


# Iterate over the different sizes of the data sets
for (n_i in 1:length(N)) {
    data_path_n_i = paste(data_path_base, N[n_i], "/X_", sep = "")
    
    # Iterate over the different realizations of the data for each size
    for (d_i in 1:10) {
        # Read in the data
        data_path = paste(data_path_n_i, d_i, ".csv", sep = "")
        X = as.matrix(read.csv(data_path, header = FALSE))
        
        ########################################################################
        # CCMM
        ########################################################################
        # Get sparse distances in dictionary of keys format
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)
        
        # Compute results CCMM
        ccmm = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE)
        timings[n_i, "CCMM"] = timings[n_i, "CCMM"] + ccmm$elapsed_time
        
        # Print the timings after CCMM
        print(timings)
        
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

        timings[n_i, "AMA"] = timings[n_i, "AMA"] +
            as.numeric(difftime(t2, t1, units = "secs"))
        
        ########################################################################
        # Store losses for n = 1000 and n = 5000
        ########################################################################
        # For n = 1000, store the values of the loss functions for comparison of
        # the accuracies of the different algorithms
        if (N[n_i] == 1000) {
            ama_cp = convert_to_cp(ama, lambdas)
            ama_loss = loss_from_cp(ama_cp$coordinates, X, w_vec, lambdas)[, 1]

            # Store the loss values in the matrix
            losses1000[, "AMA"] = losses1000[, "AMA"] + ama_loss
            losses1000[, "CCMM"] = losses1000[, "CCMM"] + ccmm$info[, "loss"]

            # Write losses to csv
            write.table(losses1000,
                        "Numerical Results/Output/n1000_losses.csv",
                        row.names = FALSE, col.names = FALSE, sep = ",")
        }
        
        # For n = 5000, store the values of the loss functions for comparison of
        # the accuracies of the different algorithms
        if (N[n_i] == 5000) {
            ama_cp = convert_to_cp(ama, lambdas)
            ama_loss = loss_from_cp(ama_cp$coordinates, X, w_vec, lambdas)
            ama_loss = ama_loss[, 1]

            # Store the loss values in the matrix
            losses5000[, "AMA"] = losses5000[, "AMA"] + ama_loss
            losses5000[, "CCMM"] = losses5000[, "CCMM"] + ccmm$info[, "loss"]

            # Write losses to csv
            write.table(losses5000,
                        "Numerical Results/Output/n5000_losses.csv",
                        row.names = FALSE, col.names = FALSE, sep = ",")
        }
        
        # Print timings to keep track of progress
        timings[n_i, "completed"] = d_i     # track data set in case of a crash
        print(timings)
        
        # Write timings to csv
        write.table(timings,
                    "Numerical Results/Output/n1000_to_5000_timings.csv",
                    row.names = FALSE, col.names = FALSE, sep = ",")
    }
}

################################################################################
# Add SSNAL timings to the table
################################################################################
timings_ssnal = read.csv(paste("Numerical Results/Output/n1000_to_5000_timings",
                               "_ssnal.csv", sep = ""), header = FALSE)
timings[, "SSNAL"] = timings_ssnal[, 1]

# Write table which now includes SSNAL timings to a csv file
write.table(timings, "Numerical Results/Output/n1000_to_5000_timings.csv",
            row.names = FALSE, col.names = TRUE, sep = ",")


################################################################################
# Compute losses for n = 1000 and n = 5000 for SSNAL
################################################################################
# Base path for the clusterpaths created by SSNAL and for the data
cp_base_path = "Numerical Results/Output/SSNAL Clusterpaths/n"
data_path_base = "Numerical Results/Data/Half Moons/Half_Moons_"

# Data set sizes
N = c(1000, 5000)

for (n_i in 1:2) {
    cp_n_i_path = paste(cp_base_path, N[n_i], "_", sep = "")
    data_path_n_i = paste(data_path_base, N[n_i], "/X_", sep = "")
    
    for (d_i in 1:10) {
        cp_path = paste(cp_n_i_path, d_i, "_ssnal_clusterpath.csv", sep = "")
        data_path = paste(data_path_n_i, d_i, ".csv", sep = "")
        
        # Load the clusterpath and the data
        ssnal_cp = as.matrix(read.csv(cp_path, header = FALSE))
        ssnal_cp = convert_to_cp(ssnal_cp, lambdas, ama = F)
        X = as.matrix(read.csv(data_path, header = FALSE))
        
        # Compute sparse weights in vectorized format
        w_vec = kernel_weights(t(X), phi)
        w_vec = knn_weights(w_vec, k, nrow(X))
        
        # Compute the loss from the clusterpaths
        ssnal_loss = loss_from_cp(ssnal_cp$coordinates, X, w_vec, lambdas)[, 1]
        
        # Add losses to appropriate columns
        if (n_i == 1) {
            losses1000[, "SSNAL"] = losses1000[, "SSNAL"] + ssnal_loss
        } else if (n_i == 2) {
            losses5000[, "SSNAL"] = losses5000[, "SSNAL"] + ssnal_loss
        }
    }
}

# Write losses to csv
write.table(losses1000,
            "Numerical Results/Output/n1000_losses.csv", 
            row.names = FALSE, col.names = TRUE, sep = ",")

write.table(losses5000,
            "Numerical Results/Output/n5000_losses.csv", 
            row.names = FALSE, col.names = TRUE, sep = ",")
