#ifndef _RSCAT_COV_H
#define _RSCAT_COV_H

#include <RcppArmadillo.h>

arma::mat calc_Sigma(const std::vector<double> alpha, const arma::mat& dist, bool usematern);
arma::mat calc_L( const arma::mat &Sigma );

arma::mat cov_powered_exponential( double sigma2, double phi, double kappa, double nugget, const arma::mat& dist);
arma::mat cov_matern( double sigma2, double phi, double nu, double nugget, const arma::mat& dist, bool uselog);
arma::mat cov_matern_vec( double sigma2, double phi, double nu, double nugget, const arma::vec& dist, bool uselog);

#endif