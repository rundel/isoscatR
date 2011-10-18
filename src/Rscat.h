#ifndef _RSCAT_H
#define _RSCAT_H

#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#include "Rscat_structs.h"


#ifdef USEMAGMA

#include <cublas.h>
#include <cuda.h>

arma::mat calc_L_gpu( std::vector<double> alpha, arma::mat& dist, bool usematern);
template <typename T> void init_gpu_mem(T *source, T *d_B, int n);
template <typename T> void clean_gpu_mem(T *dest, T *d_B, int n);

void mag_schol_debug(arma::mat &B, bool gpu = true);
RcppExport SEXP magma_chol(SEXP rX, SEXP rGPU, SEXP rFLOAT);

void checkCudaError(const char *msg);
void checkCublasError(const char *msg);
std::string cublasGetErrorString(cublasStatus err);

#endif


// Locate functions
void update_Location(GlobalParams &p, GlobalOptions &opt);
void update_LocationCV(GlobalParams &p, GlobalOptions &opt);
void init_locate(GlobalParams &p, GlobalOptions &opt);



// init functions
void parseArgs(GlobalParams &p, SEXP rBoundary, SEXP rLocations, SEXP rRegionNames,
                                SEXP rGenotypes, SEXP rIndivID, SEXP rNalleles);
void init_params(GlobalParams &p, GlobalOptions &opt);
void init_proposal_sd(GlobalParams &p, GlobalOptions &opt);
void init_attempts(GlobalParams &p);
void calc_counts(GlobalParams &p, GlobalOptions &opt);
void calc_params(GlobalParams &p, GlobalOptions &opt);
void open_allelefiles(GlobalParams &p, GlobalOptions &opt);
void close_allelefiles(GlobalParams &p, GlobalOptions &opt);
void open_cvfiles(GlobalParams &p, GlobalOptions &opt);
void close_cvfiles(GlobalParams &p, GlobalOptions &opt);

//MCMC Functions

double tuneScale(double acc_rate);
void MCMCTune(GlobalParams &p, GlobalOptions &opt);

Rcpp::List MCMCLoop(GlobalParams &p, GlobalOptions &opt, int Niter, int Nthin, bool burnin, bool thin);
void MCMCStep( GlobalParams &p, GlobalOptions &opt, bool burnin);    
void update_Beta(GlobalParams &p, GlobalOptions &opt);
void update_X(GlobalParams &p, GlobalOptions &opt);
void update_Alpha(GlobalParams &p, GlobalOptions &opt);
void update_Mu(GlobalParams &p, GlobalOptions &opt);
void update_Eta(GlobalParams &p, GlobalOptions &opt);
void update_Xi(GlobalParams &p, GlobalOptions &opt);
void update_anisotropy(GlobalParams &p, GlobalOptions &opt);


// Option Functions
void parseOptions(SEXP sexpRopt, GlobalOptions &opt);


// Output Functions
void outputAccepts(GlobalParams &p, GlobalOptions &opt);
void outputTuning(GlobalParams &p, GlobalOptions &opt);


// cov Functions

arma::mat calc_L( std::vector<double> alpha, arma::mat& dist, bool usematern);
arma::mat cov_powered_exponential( double sigma2, double phi, double kappa, double nugget, arma::mat& dist);
arma::mat cov_matern( double sigma2, double phi, double nu, double nugget, arma::mat& dist, bool uselog);
arma::mat cov_matern_vec( double sigma2, double phi, double nu, double nugget, arma::vec& dist, bool uselog);
RcppExport SEXP R_cov_matern( SEXP rsigma2, SEXP rphi, SEXP rnu, SEXP rnugget, SEXP rdist, SEXP distmat, SEXP ruselog);



// utility functions
RcppExport SEXP read_allele_file(SEXP Rfile);
RcppExport SEXP R_calc_distance(SEXP x, SEXP y);
RcppExport SEXP R_calc_distance_to_point(SEXP px, SEXP py, SEXP x, SEXP y);
RcppExport SEXP prec_sum(SEXP Rvec);

double calc_accept_ratio(unsigned int accept, unsigned int attempt);

arma::mat calc_rotation_mat(double angle);
arma::mat calc_stretch_mat(double ratio);

arma::colvec calc_LogLik(arma::mat theta, arma::colvec sumExpTheta, arma::mat count, arma::colvec sumCount);
arma::colvec calc_multinom_loglik( const arma::mat theta, const arma::mat count, arma::colvec sumCount);

arma::mat calc_theta(const double &mu, const arma::rowvec &eta, const double &xi, const arma::mat &L, const arma::mat &X);
arma::mat calc_f(const arma::mat &theta);

arma::mat calc_distance_mat(arma::mat const &m);
arma::mat calc_distance_mat(arma::colvec const &xc, arma::colvec const &yc);

double calc_distance(double x1, double y1, double x2, double y2);

arma::colvec dnorm(arma::colvec x);
double dnorm(double x);

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2);
int isInsideBoundary( double x, double y, arma::mat boundary);


// main functions

RcppExport SEXP mcmc_main(SEXP rChain,
              SEXP rBoundary,   // px2 matrix
              SEXP rLocations,  // Rx2 matrix
              SEXP rGenotypes, // 2Ix(L+1) matrix
              SEXP rIndRegion,
              SEXP rNalleles,
              SEXP rNiter,
              SEXP rNthin,
              SEXP rNburn,
              SEXP rCVIndivs,
              SEXP rCVGenotypes,
              SEXP rOpt );


#endif
