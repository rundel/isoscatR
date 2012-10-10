#include <Rmath.h>
#include <iostream>

#include "scatR_cov.h"

arma::mat calc_Sigma(const std::vector<double> alpha, const arma::mat& dist, bool usematern)
{
    return (usematern) 
            ? cov_matern(alpha[0], alpha[1], alpha[2], alpha[3], dist, false)
            : cov_powered_exponential(alpha[0], alpha[1], alpha[2], alpha[3], dist);
    
}

arma::mat calc_L(const arma::mat& Sigma)
{
    arma::mat L;
    try {
        L = arma::chol(Sigma).t(); // Chol returns Upper Triangular, return Low Tri
    } catch(...) {
        L.reset();
    }
    
    if (L.is_empty()) {
        std::cout << "Cholesky decomposition failed! ";
    }

    return L;
}


// Powered exponential covariance function with sill and nugget effect
arma::mat cov_powered_exponential( double sigma2, double phi, double kappa, double nugget, const arma::mat& dist)
{
    return sigma2 * arma::exp( -arma::pow(dist / phi, kappa) ) 
           + arma::eye<arma::mat>(dist.n_rows,dist.n_cols)*nugget;
}

// Matern covariance function with sill and nugget effect
arma::mat cov_matern(double sigma2, double phi, double nu, double nugget, const arma::mat& dist, bool uselog)
{
    int nr = dist.n_rows;
    
    double diag = (uselog) ? log(sigma2+nugget) : sigma2+nugget;
    arma::mat S = arma::eye<arma::mat>(nr,nr) * diag;
    
    for(int i=0; i<nr-1; i++) {
        for(int j=i+1; j<nr; j++) {
            
            if (dist(i,j) > 600 * phi) {
                S(i,j) = 0;
                S(j,i) = 0;
                continue;
            }
            
            double temp = dist(i,j) / phi; 
            
            S(i,j) = (uselog)
                     ? (log(sigma2) - Rf_lgammafn(nu) - (nu-1)*log(2) + nu*log(temp)+log(Rf_bessel_k(temp,nu,1)))
                     : (sigma2 / (Rf_gammafn(nu) * pow(2,(nu-1)) ) * pow(temp,nu) * Rf_bessel_k(temp,nu,1));
            S(j,i) = S(i,j);
        }
    }
    
    return( (uselog) ? exp(S) : S );
}

arma::mat cov_matern_vec( double sigma2, double phi, double nu, double nugget, const arma::vec& dist, bool uselog)
{
    arma::vec cov(dist.n_elem);
    
    for(int i=0; i<dist.n_elem; i++) {
        if (dist(i) == 0) {
            cov(i) = (uselog) ? log(sigma2+nugget) : sigma2+nugget;
        } else if (dist(i) > 600 * phi) {
            cov(i) = 0;
        } else {
            double temp = dist(i) / phi; 
            
            cov(i) = (uselog) ? (log(sigma2) - Rf_lgammafn(nu) - (nu-1)*log(2) + nu*log(temp)+log(Rf_bessel_k(temp,nu,1)))
                              : (sigma2 / (Rf_gammafn(nu) * pow(2,(nu-1)) ) * pow(temp,nu) * Rf_bessel_k(temp,nu,1));    
        }
    }
    
    return( (uselog) ? exp(cov) : cov );
}
