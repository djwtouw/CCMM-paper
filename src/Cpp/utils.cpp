//[[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include <RcppArmadillo.h>


double loss(arma::mat& A, arma::vec& w, arma::mat& X, double lambda)
{
    double result = 0.5*pow(arma::norm(X-A, "fro"), 2);
    int iterator = 0;
    int n = A.n_rows;
    
    for (int i = 0; i < (n - 1); i++) {
        for (int j = i + 1; j < n; j++) {
            result += lambda * w(iterator) * arma::norm(A.row(i) - A.row(j), 2);
            iterator++;
        }
    }
    
    return result;
}


// [[Rcpp::export]]
arma::vec loss_from_cp(arma::mat& results, arma::mat& X, arma::vec& w, arma::vec& lambdas)
{
    int n = X.n_rows;
    int p = X.n_cols;
    int n_lambdas = lambdas.n_elem;
    arma::mat A(X);
    arma::vec losses(n_lambdas);
    
    for (int l = 0; l < n_lambdas; l++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                A(i, j) = results(l * n + i, j);
            }
        }
        
        losses(l) = loss(A, w, X, lambdas(l));
    }
    
    return losses;
}
