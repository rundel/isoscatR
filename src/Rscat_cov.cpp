#include <Rmath.h>
#include "Rscat.h"


using namespace arma;
using namespace Rcpp;
using namespace std;

mat calc_L( vector<double> alpha, mat& dist, bool usematern) {
    
    mat m = (usematern) 
            ? cov_matern(alpha[0], alpha[1], alpha[2], alpha[3], dist, false)
            : cov_powered_exponential(alpha[0], alpha[1], alpha[2], alpha[3], dist);
    
    mat L = trans( chol(m) ); // Chol returns Upper Triangular, return LT

    if (L.is_empty()) {
        cout << "Cholesky decomposition failed! ";
        cout << "Alpha: " << alpha[0] << " " << alpha[1] << " "
                          << alpha[2] << " " << alpha[3] << endl;
    }

    return(L);
}


// Powered exponential covariance function with sill and nugget effect
mat cov_powered_exponential( double sigma2, double phi, double kappa, double nugget, mat& dist) {
    
    int nr = dist.n_rows;
    mat L = sigma2 * exp( -pow(dist / phi, kappa) ) + eye<mat>(nr,nr)*nugget;
    
    return( L );
}

// Matern covariance function with sill and nugget effect
mat cov_matern( double sigma2, double phi, double nu, double nugget, mat& dist, bool uselog) {
    
    int nr = dist.n_rows;
    
    double diag = (uselog) ? log(sigma2+nugget) : sigma2+nugget;
    mat L = eye<mat>(nr,nr) * diag;
    
    for(int i=0; i<nr-1; i++) {

        for(int j=i+1; j<nr; j++) {
            
            if (dist(i,j) > 600 * phi) {
                L(i,j) = 0;
                L(j,i) = 0;
                continue;
            }
            
            double temp = dist(i,j) / phi; 
            
            L(i,j) = (uselog)
                     ? (log(sigma2) - Rf_lgammafn(nu) - (nu-1)*log(2) + nu*log(temp)+log(Rf_bessel_k(temp,nu,1)))
                     : (sigma2 / (Rf_gammafn(nu) * pow(2,(nu-1)) ) * pow(temp,nu) * Rf_bessel_k(temp,nu,1));
            L(j,i) = L(i,j);
        }
    }
    
    return( (uselog) ? exp(L) : L );
}

mat cov_matern_vec( double sigma2, double phi, double nu, double nugget, vec& dist, bool uselog) {
    
    vec cov(dist.n_elem);
    
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


SEXP R_cov_matern( SEXP rsigma2, SEXP rphi, SEXP rnu, SEXP rnugget, SEXP rdist, SEXP rdistmat, SEXP ruselog) {
     
    double nugget = as<double>(rnugget);
    double sigma2 = as<double>(rsigma2);
    double nu     = as<double>(rnu);
    double phi    = as<double>(rphi);
    bool distmat  = as<bool>(rdistmat);
    bool uselog   = as<bool>(ruselog);
    
    SEXP res;
    if (distmat) {
        mat dist = as<mat>(rdist);
        
        res = wrap(cov_matern(sigma2, nu, phi, nugget, dist, uselog));
    } else {
        //NumericVector tDist(rdist);
        vec dist = as<vec>(rdist);
        
        res = wrap(cov_matern_vec(sigma2, nu, phi, nugget, dist, uselog));
    }
    
    return( res );
}