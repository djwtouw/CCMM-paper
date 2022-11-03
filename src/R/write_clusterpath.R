write_clusterpath <- function(clusterpath, path)
{
    data = cbind(1, cbind(1, clusterpath$coordinates))
    data[, 1] = rep(c(1:clusterpath$n), times = length(lambdas))
    data[, 2] = rep(lambdas, each = clusterpath$n)
    
    write.table(data, path, row.names = FALSE, col.names = FALSE)
}