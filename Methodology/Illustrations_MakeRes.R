################################################################################
# Compute the results for dense versus sparse weights comparison and the effect
# of the symmetric circulant on the clusterpaths.
################################################################################

rm(list = ls())
library(CCMMR)
source("src/R/write_clusterpath.R")


################################################################################
# Dense vs Sparse
################################################################################

X = as.matrix(read.csv("Methodology/Data/Dense_vs_Sparse/2half_moons_X.csv", 
                       header = FALSE))
y = as.matrix(read.csv("Methodology/Data/Dense_vs_Sparse/2half_moons_y.csv", 
                       header = FALSE))

### UNIT WEIGHTS
# Get dense unit weights
W = sparse_weights(X, nrow(X) - 1, 0.0)

# Sequence for lambda
lambdas = seq(0, 20, 0.05)

# Results
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Dense_vs_Sparse/res_unit.csv")


### DENSE WEIGHTS
# Get dense weights in dictionary of keys format
W = sparse_weights(X, nrow(X) - 1, 4.0)

# Sequence for lambda
lambdas = seq(0, 70, 0.25)

# Results
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Dense_vs_Sparse/res_dense.csv")


### SPARSE WEIGHTS
# Get sparse weights in dictionary of keys format with k = 10
W = sparse_weights(X, 10, 4.0)

# Set a sequence for lambda
lambdas = seq(0, 240, 0.5)

# Perform clustering
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Dense_vs_Sparse/res_sparse.csv")


################################################################################
# Connectedness
################################################################################

X = as.matrix(read.csv("Methodology/Data/Connectedness/X.csv", 
                       header = FALSE))
y = as.matrix(read.csv("Methodology/Data/Connectedness/y.csv", 
                       header = FALSE))


# WITH SYMMETRIC CIRCULANT
# Compute sparse weights
W = sparse_weights(X, 3, 1.0, connected = TRUE, connection_type = "SC")

# Sequence for lambda
lambdas = seq(0, 130, 0.5)

# Compute results
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Connectedness/res_SC.csv")
write.table(cbind(W$keys, W$values), 
            "Methodology/Output/Connectedness/knn_weights_SC.csv",
            row.names = FALSE, col.names = FALSE)

# WITH MINIMUM SPANNING TREE
# Compute sparse weights
W = sparse_weights(X, 3, 1.0, connected = TRUE, connection_type = "MST")

# Sequence for lambda
lambdas = seq(0, 1030, 0.5)

# Compute results
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Connectedness/res_MST.csv")
write.table(cbind(W$keys, W$values), 
            "Methodology/Output/Connectedness/knn_weights_MST.csv",
            row.names = FALSE, col.names = FALSE)

# DISCONNECTED
# Compute sparse weights
W = sparse_weights(X, 3, 1.0, connected = FALSE)

# Sequence for lambda
lambdas = seq(0, 160, 0.5)

# Compute results
res = convex_clusterpath(X, W, lambdas)
plot(res, y + 1)

write_clusterpath(res, "Methodology/Output/Connectedness/res_DC.csv")
write.table(cbind(W$keys, W$values),
            "Methodology/Output/Connectedness/knn_weights_DC.csv",
            row.names = FALSE, col.names = FALSE)


################################################################################
# Dendrogram
################################################################################

# Generate data
set.seed(6)
X = matrix(rnorm(14), ncol = 2)
y = rep(1, nrow(X))

# Get sparse distances in dictionary of keys format with k = 3
W = sparse_weights(X, 3, 4.0, connected = FALSE)

# Sequence for lambda
lambdas = seq(0, 45, 0.01)

# Compute results
res = convex_clusterpath(X, W, lambdas)
plot(res, y)

# Generate hclust object
hcl = as.hclust(res)
hcl$height = sqrt(hcl$height)

# Plot clusterpath and dendrogram
par(mfrow=c(1,2))
plot(res, y, label = c(1:nrow(X)))
plot(hcl, ylab = expression(sqrt(lambda)), xlab = NA, sub = NA, main = NA)

# Write results to csv
write_clusterpath(res, "Methodology/Output/Dendrogram/res.csv")
write.table(hcl$merge, "Methodology/Output/Dendrogram/merge.csv",
            row.names = FALSE, col.names = FALSE)
write.table(hcl$height, "Methodology/Output/Dendrogram/height.csv",
            row.names = FALSE, col.names = FALSE)
write.table(hcl$order, "Methodology/Output/Dendrogram/order.csv",
            row.names = FALSE, col.names = FALSE)
