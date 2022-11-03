convert_to_cp <- function(ama_or_ssnal, lambdas, ama = TRUE) {
  n_lambdas = length(lambdas)
  
  if (ama) {
    n = ncol(ama_or_ssnal$U[[1]])
    p = nrow(ama_or_ssnal$U[[1]])
    
    coordinates = (matrix(NA, nrow = n * n_lambdas, ncol = p))
    
    for (i in 1:n_lambdas) {
      temp = t(ama_or_ssnal$U[[i]])
      coordinates[((i - 1) * n + 1):(i * n), ] = temp
    }
    
  } else {
    n = nrow(ama_or_ssnal) / n_lambdas
    p = ncol(ama_or_ssnal) - 2
    coordinates = ama_or_ssnal[, 3:(p + 2)]
  }
  
  result = list()
  result$coordinates = coordinates
  result$lambdas = lambdas
  result$n = n
  class(result) = "cvxclust"
  
  return(result)
}