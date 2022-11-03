################################################################################
# Script that verifies if the output of the different algorithms agree given 
# that the inputs to all are the same.
################################################################################

rm(list = ls())

library(CCMMR)
library(cvxclustr)

source("src/R/convert_to_cp.R")
Rcpp::sourceCpp("src/Cpp/utils.cpp")


# Load data
X = as.matrix(read.csv(paste("Numerical Results/Data/Half Moons/",
                             "Half_Moons_0200/X.csv", sep = ""), 
                       header = FALSE))
y = rep(1, nrow(X))

# Set a sequence for lambda
lambdas = seq(0, 10, 0.2)

################################################################################
# Run the CCMM algorithms for n = 200
################################################################################

# Get sparse distances in dictionary of keys format with k = 10
W = sparse_weights(X, 10, 2.0, scale = FALSE, connected = FALSE)

# Compute results CCMM
ccmm = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE)
plot(ccmm, y)


################################################################################
# Run AMA for n = 200
################################################################################

w_vec = kernel_weights(t(X), 2)
w_vec = knn_weights(w_vec, 10, nrow(X))

nu = AMA_step_size(w_vec, nrow(X))
ama = cvxclust_path_ama(t(X), w_vec, lambdas, nu = nu, tol = 1e-3)
ama = convert_to_cp(ama, lambdas)

plot(ama, y)


################################################################################
# Load the results of SSNAL for n = 200
################################################################################

ssnal = as.matrix(read.csv(paste("Numerical Results/Output/SSNAL Clusterpaths/", 
                                 "n200_ssnal_clusterpath.csv",
                                 sep = ""), header = FALSE))
ssnal = convert_to_cp(ssnal, lambdas, ama = F)
plot(ssnal, y)


################################################################################
# After visual inspection of the clusterpaths, also check if the values of the
# obtained loss function is similar for each algorithm
################################################################################

# Select the loss computed by each of the algorithms, ignoring the loss for 
# lambda = 0
loss_ccmm = ccmm$info$loss[-1]
loss_ama = loss_from_cp(ama$coordinates, X, w_vec, lambdas)[-1, 1]
loss_ssnal = loss_from_cp(ssnal$coordinates, X, w_vec, lambdas)[-1, 1]

# Plot the relative loss with the values obtained by SSNAL as benchmark
plot(lambdas[-1], loss_ccmm / loss_ssnal, type = "l", col = "blue", 
     ylim = c(0.99999, 1.00015), ylab = "Relative Loss w.r.t. SSNAL")
lines(lambdas[-1], loss_ama / loss_ssnal, col = "red")
lines(lambdas[-1], loss_ssnal / loss_ssnal, col = "black")
legend("topleft", legend=c("CCMM", "AMA", "SSNAL"),
       col=c("blue", "red", "black"), lty = 1, cex=0.8)
