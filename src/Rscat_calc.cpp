#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

mat calc_rotation_mat(double angle) {
    
    mat res(2,2);
    res(0,0) =  cos(angle);
    res(0,1) = -sin(angle);
    res(1,0) =  sin(angle);
    res(1,1) =  cos(angle);
    
    return(res);
}

mat calc_stretch_mat(double ratio) {
    
    mat res(2,2);
    res.eye();
    res(1,1) = 1/ratio;
    
    return(res);
}

double calc_multinomial_loglik(const mat &theta, const mat &count, const colvec &sumCount) {
    
    return accu( sum(count % theta,1) - sumCount % log(sum(exp(theta),1)) );
}

double calc_multivar_normal_loglik(const mat &x, const mat &mu, const mat &sigmaInv, double sigmaDet) {

    double res = 0.0;
    for(int j=0; j<x.n_cols; j++) {
        mat tmp = -0.5*( log(sigmaDet) + trans(x.col(j) - mu.col(j)) * sigmaInv * (x.col(j) - mu.col(j)) );
        //tmp.print();
        res += accu(tmp);
    }
        
    return res;
}

mat calc_f(const mat &theta) {
    
    int nAlleles = theta.n_cols;
    
    mat expTheta = exp(theta);
    mat sumExpTheta = sum(expTheta,1) * ones<rowvec>(nAlleles);
    
    return( expTheta / sumExpTheta );
}