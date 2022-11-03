################################################################################
# Script for analyzing the gamma data set.
################################################################################

# Clear environment
rm(list = ls())

# Load required packages
library(dslabs)
library(uwot)
library(CCMMR)
library(viridis)
library(cvxclustr)

source("src/R/min_clusters.R")

# Reset plot viewer
par(mfrow = c(1, 1))


# Check if there is a file with the results, if not, create one, if so, make
# sure the results for the current data set can be added
current_data = "gamma"
if (!(file.exists("Numerical Results/Output/uci_results.csv"))) {
    results = matrix(NA, nrow = 1, ncol = 14)
    colnames(results) = c("n", "p", "c", "n_min", "n_max", "p_used", "k", 
                          "phi", "max_lambda", "n_lambdas", "time_ccmm", 
                          "time_ama", "time_ssnal", "min_clusters")
    rownames(results) = current_data
} else {
    results = as.matrix(read.csv("Numerical Results/Output/uci_results.csv"))
    
    if (!(current_data %in% rownames(results))) {
        results = rbind(results, NA)
        rownames(results)[nrow(results)] = current_data
    }
}


# Load and preprocess the data
data = read.csv(paste("Numerical Results/Data/UCI/", current_data, ".csv", 
                      sep = ""), sep = ",", header = FALSE)

# Separate X and y
X = as.matrix(data[, -11])
y = data[, 11]
y = as.numeric(as.factor(y))
y = y - min(y) + 1

# Remove unnecessary variables
rm(data)

# If it exists, load preprocessed data
if (!file.exists(paste("Numerical Results/Data/UCI/", current_data,
                       "_X_umap.csv", sep = ""))) {
    
    # Apply shuffling
    set.seed(1)
    indices = sample(1:nrow(X))
    X = X[indices, ]
    y = y[indices]
    
    # Apply preprocessing
    X.pre <- umap(X, n_neighbors = 50, min_dist = 0.001, verbose = TRUE,
                  n_components = length(unique(y)))
    
    # Ensure the scale of the data is not too large, that can cause nonzero
    # weights that are numerically indistinguishable from zero
    X.pre = X.pre / (0.5 * sd(X.pre))
    
    # Store X
    write.table(X.pre, paste("Numerical Results/Data/UCI/", current_data, 
                             "_X_umap.csv", sep = ""),
                sep = ",", col.names = FALSE, row.names = FALSE)
    
    # Store y
    write.table(y, paste("Numerical Results/Data/UCI/", current_data, 
                         "_y_umap.csv", sep = ""),
                sep = ",", col.names = FALSE, row.names = FALSE)
    
    rm(indices)
} else {
    # Load preprocessed data
    X.pre = read.csv(paste("Numerical Results/Data/UCI/", current_data, 
                           "_X_umap.csv", sep = ""), sep = ",", header = FALSE)
    X.pre = as.matrix(X.pre)
    y = read.csv(paste("Numerical Results/Data/UCI/", current_data, 
                       "_y_umap.csv", sep = ""), sep = ",", header = FALSE)
    y = as.matrix(y)[, 1]
}

# Plot the transformed data
plot(X.pre[, 1:2], pch = 19, col = cividis(max(y))[y])

# Use AMA if system allows it
use_ama = FALSE # <- This data set is too large

# Settings for the weights
k = 50
phi = 2

# Set W and target range of clusters
W = sparse_weights(X.pre, k = k, phi = phi, scale = FALSE, connected = FALSE)
min(W$values)

# Minimum number of clusters that is possible to get with this weight matrix
min_c = min_clusters(W$keys)

# Find the value for lambda that corresponds to the minimum number of clusters
res = convex_clustering(X.pre, W, target_low = min_c, target_high = 50, 
                        scale = FALSE, center = FALSE)

# Create a sequence of lambdas used for the clusterpath
max_lambda = res$info[, "lambda"][nrow(res$info)] * 1.05
lambdas = max_lambda / 199^3 * c(0:199)^3

# Store lambdas
write.table(lambdas, paste("Numerical Results/Data/UCI/", current_data, 
                           "_lambdas.csv", sep = ""),
            sep = ",", col.names = FALSE, row.names = FALSE)

# Compute the clusterpath using CCMM
res = convex_clusterpath(X.pre, W, lambdas, scale = FALSE, center = FALSE)

# Plot the clusterpath
if (ncol(X.pre) == 2) {
    par(mfrow = c(1, 1))
    plot(res, col = cividis(max(y))[y])
}

# Fill results
results[current_data, "n"] = nrow(X.pre)
results[current_data, "p"] = ncol(X)
results[current_data, "c"] = length(unique(y))
results[current_data, "n_min"] = min(table(y))
results[current_data, "n_max"] = max(table(y))
results[current_data, "p_used"] = ncol(X.pre)
results[current_data, "k"] = k
results[current_data, "phi"] = phi
results[current_data, "time_ccmm"] = res$elapsed_time
results[current_data, "max_lambda"] = max_lambda
results[current_data, "n_lambdas"] = length(lambdas)
results[current_data, "min_clusters"] = res$info[length(lambdas), "clusters"]

# Write results to file
write.table(results, "Numerical Results/Output/uci_results.csv",
            col.names = TRUE, row.names = TRUE, sep = ",")
print(results)

# If possible, also use AMA
if (use_ama) {
    # Get the weights in vector format
    w_vec = kernel_weights(t(X.pre), phi)
    w_vec = knn_weights(w_vec, k, nrow(X.pre))
    
    # Compute recommended step size nu
    nu = AMA_step_size(w_vec, nrow(X.pre))
    
    # Compute AMA result
    t1 = Sys.time()
    ama = cvxclust_path_ama(t(X.pre), w_vec, lambdas, nu = nu)
    t2 = Sys.time()
    
    results[current_data, "time_ama"] = 
        as.numeric(difftime(t2, t1, units = "secs"))
}

# Write results to file
write.table(results, "Numerical Results/Output/uci_results.csv",
            col.names = TRUE, row.names = TRUE, sep = ",")
print(results)
