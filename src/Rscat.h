#ifndef _Rscat_RCPP_HELLO_WORLD_H
#define _Rscat_RCPP_HELLO_WORLD_H

#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


struct GlobalParams {
    
    int chain_num;
    
    int nRegions;
    int nLoci;
    int nInd;
    
    Rcpp::IntegerVector         nAlleles;
    
    arma::imat                  genotypes;
    Rcpp::IntegerVector         cvIndivs;
    arma::imat                  cvGeneotypes;
    
    Rcpp::IntegerVector         indRegion;
    
    arma::mat                   predDist;
    
    std::vector<arma::mat>      count;          // NA  RxA[l]
    std::vector<arma::colvec>   sumCount;       // NA  Rx1
    std::vector<arma::mat>      allele_freq;    // NA  RxA[l]
    
    std::vector<double>         mu;             // NA
    std::vector<arma::rowvec>   eta;            // NA  1xA[l]
    std::vector<double>         beta;           // 1
    std::vector<double>         xi;             // 1
    
    std::vector<double>         alpha;          // 4
    
    std::vector<arma::mat>      X;              // NA  RxA[l]
    std::vector<arma::mat>      theta;          // NA  RxA[l]
    
    double                      anisoRatio;
    double                      anisoAngle;
    
    arma::mat                   anisoRatioMat;
    arma::mat                   anisoAngleMat;
    
    
    arma::mat                   dist;           // RxR
    arma::mat                   L;              // RxR
    
    std::vector<arma::colvec>   logLik;         // NA  Rx1

    arma::mat                   locs;           // Rx2
    arma::mat                   locsTrans;      // Rx2

    
    arma::mat                   boundary;
    
    std::vector<arma::colvec>   x_sd;        // NA   Rx1
    arma::rowvec                mu_sd;       
    arma::rowvec                eta_sd;      
    arma::rowvec                xi_sd;
    arma::rowvec                beta_sd;
    arma::rowvec                alpha_sd;
    double                      angle_sd;
    double                      ratio_sd;
    
    std::vector<arma::ucolvec>   Xattempt;        // NA   Rx1
    std::vector<arma::ucolvec>   Xaccept;         // NA   Rx1
    
    arma::urowvec                muAttempt;       // NA   
    arma::urowvec                muAccept;        // NA 
    
    arma::urowvec                etaAttempt;       
    arma::urowvec                etaAccept;       
    
    arma::urowvec                xiAttempt;
    arma::urowvec                xiAccept;
    
    arma::urowvec                betaAttempt;
    arma::urowvec                betaAccept;
    
    arma::urowvec                alphaAttempt;
    arma::urowvec                alphaAccept;
    
    unsigned int                 angleAttempt;
    unsigned int                 angleAccept;
    
    unsigned int                 ratioAttempt;
    unsigned int                 ratioAccept;

    std::vector<std::vector<boost::iostreams::filtering_ostream*> > alfileGzStreams;
    std::vector<std::vector<std::ofstream*> > alfileStreams;
    
    std::vector<boost::iostreams::filtering_ostream*> cvfileGzStreams;
    std::vector<std::ofstream*> cvfileStreams;
};

struct GlobalOptions {
    bool VERBOSE;
    std::string TMPDIR;
    std::string FILESUFFIX;
    
    bool ADAPT;
    int TUNEINTERVAL;
    
    bool USEMATERN;
    
    bool RETURNFIT;
    
    bool LOCATE;
    double MAXCELL;
    
    bool CROSSVALIDATE;
    
    double PSEUDOCOUNT;
    
    Rcpp::LogicalVector FIXALPHA;
    Rcpp::NumericVector ALPHA;
    Rcpp::NumericVector ALPHAMIN;
    Rcpp::NumericVector ALPHAMAX;
    
    Rcpp::NumericVector ALPHASD;
    
    bool FIXANGLE;
    double ANGLE;
    double ANGLESD;
    
    bool FIXRATIO;
    double RATIO;
    double RATIOSD;
    
    
    Rcpp::NumericVector XIRANGE;
    bool FIXXI;
    Rcpp::NumericVector XI;
    double XISD;
    double SIGMAXI;
    
    bool FIXMU;
    Rcpp::NumericVector MU;
    double MUSD;
    
    bool FIXETA;
    Rcpp::NumericVector ETA;
    double ETASD;
    
    Rcpp::NumericVector BETARANGE;
    bool FIXBETA;
    Rcpp::NumericVector BETA;
    double BETASD;
    double SIGMABETA;

    double XSD;
    
    double LOCALSD;
    double GLOBALSD;
    
    double NULLPROB;
    double DELTA;
};

// Locate functions
void update_Location(GlobalParams &p, GlobalOptions &opt);
void update_LocationCV(GlobalParams &p, GlobalOptions &opt);
void init_locate(GlobalParams &p, GlobalOptions &opt);



// init functions
void parseArgs(GlobalParams &p, SEXP rBoundary, SEXP rLocations, SEXP rRegionNames,
                                SEXP rGeneotypes, SEXP rIndivID, SEXP rNalleles);
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

arma::mat calc_L( std::vector<double> alpha, arma::mat dist, bool usematern);
arma::mat cov_powered_exponential( double sigma2, double phi, double kappa, double nugget, arma::mat dist);
arma::mat cov_matern( double sigma2, double phi, double nu, double nugget, arma::mat dist, bool uselog);
arma::mat cov_matern_vec( double sigma2, double phi, double nu, double nugget, arma::vec dist, bool uselog);
RcppExport SEXP R_cov_matern( SEXP rsigma2, SEXP rphi, SEXP rnu, SEXP rnugget, SEXP rdist, SEXP distmat, SEXP ruselog);



// utility functions
RcppExport SEXP read_allelefile(SEXP Rfile);
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
Rcpp::NumericVector calc_distance_to_point(double px, double py, arma::colvec &xc, arma::colvec &yc);

arma::colvec dnorm(arma::colvec x);
double dnorm(double x);

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2);
int isInsideBoundary( double x, double y, arma::mat boundary);


// main functions

RcppExport SEXP mcmc_main(SEXP rChain,
              SEXP rBoundary,   // px2 matrix
              SEXP rLocations,  // Rx2 matrix
              SEXP rGeneotypes, // 2Ix(L+1) matrix
              SEXP rIndRegion,
              SEXP rNalleles,
              SEXP rNiter,
              SEXP rNthin,
              SEXP rNburn,
              SEXP rCVIndivs,
              SEXP rCVGeneotypes,
              SEXP rOpt );


#endif
