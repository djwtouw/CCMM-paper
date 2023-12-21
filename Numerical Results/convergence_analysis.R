################################################################################
# Script to analyze the convergence of CCMM. First part estimates the order of
# convergence. Second part plots loss function value progression for sequences
# that were long enough.
################################################################################

rm(list = ls())

# Load libraries
library(CCMMR)

# Source helper functions
source("src/R/save_load.R")


# Set parameters for the weights
k = 15
phi = 2
lambdas = seq(0, 110, 0.2)

# List for all loss values obtained during optimization
losses = list()
convergence_norms = list()

# The base path to the folder where the data is stored
data_path_base = "Numerical Results/Data/Half Moons/Half_Moons_"

for (n_i in c("1000", "2000", "3000", "4000", "5000")){
    data_path_n_i = paste(data_path_base, n_i, "/X_", sep = "")

    # Iterate over the different realizations of the data for each size
    for (d_i in 1:10) {
        cat("n = ", n_i, ", r = ", d_i, "\n", sep = "")
        # Read in the data
        data_path = paste(data_path_n_i, d_i, ".csv", sep = "")
        X = as.matrix(read.csv(data_path, header = FALSE))
        
        # Get sparse distances in dictionary of keys format
        W = sparse_weights(X, k = k, phi = phi, scale = FALSE, 
                           connected = FALSE)
        
        # Compute results CCMM
        ccmm = convex_clusterpath(X, W, lambdas, scale = FALSE, center = FALSE, 
                                  save_losses = TRUE,
                                  save_convergence_norms = TRUE)
        
        # Add losses
        losses[[length(losses) + 1]] = ccmm$losses
        
        # Add convergence norms
        convergence_norms[[length(convergence_norms) + 1]] =
            ccmm$convergence_norms
    }
}

# Estimates for the order of convergence
orders = list()

for (convergence_norms_i in convergence_norms) {
    for (x in convergence_norms_i) {
        if (length(x) < 3) {
            next
        }
        
        # Container for orders
        orders_i = rep(0, length(x) - 2)
        
        # Compute alphas
        for (n in 3:length(x)) {
            N = log(x[n] / x[n - 1])
            D = log(x[n - 1] / x[n - 2])
            orders_i[n - 2] = N / D
        }
        
        # Remove -Infs due to consecutive iterates being equal (indicating
        # actual convergence)
        orders_i = orders_i[orders_i > -Inf]
        
        # Add
        orders[[length(orders) + 1]] = orders_i
    }
}

# Create histogram
x = unlist(orders)
hist(x[((x > 0) * (x < 2)) == 1], breaks = 1000)
mean(x)
median(x)

# Add vertical lines for median and mean
abline(v = median(x), col = "red", lty = 2, lwd = 2)
abline(v = mean(x), col = "blue", lty = 3, lwd = 2) 

# Percentage of the data featured in the histogram
round(length(x[((x > 0) * (x < 2)) == 1]) / length(x) * 100, 1)

# Percentage of the estimates below one
round(length(x[x < 1]) / length(x) * 100, 1)

# Store results
save_bin(x, "Numerical Results/Output/convergence_analysis.bin")


## Analysis of loss function value progression for minimizations that took at
## least 100 iterations to converge
# Get losses for all minimizations that took 100 iterations
losses_100 = list()
for (losses_i in losses) {
    for (x in losses_i) {
        if (length(x) < 101) {
            next
        }
        
        x = x[c(1:101)]
        x = (x - min(x)) / (max(x) - min(x))
        
        losses_100[[length(losses_100) + 1]] = x
    }
}

# Plot all vectors in the same plot with alpha = 0.5
plot(NULL, xlim = c(1, max(sapply(losses_100, length))), 
     ylim = c(min(unlist(losses_100)), max(unlist(losses_100))), type = "n")

# In the mean time, also compute the mean values of the normalized losses
mean_loss = rep(0, 101)
for (x in losses_100) {
    lines(x, col = rgb(1, 0, 0, 0.05))
    mean_loss = mean_loss + x
}
mean_loss = mean_loss / length(losses_100)
lines(mean_loss)

# Write losses to disc
losses_100 = matrix(unlist(losses_100), nrow = 101, byrow = FALSE)
write.table(
    losses_100, "Numerical Results/Output/convergence_analysis_losses.csv",
    sep = ",", row.names = FALSE, col.names = FALSE
)
