#ifndef _RSCAT_STRUCTS_H
#define _RSCAT_STRUCTS_H

#include <RcppArmadillo.h>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#ifdef USEMAGMA
#include "Rscat_gpu.h"
#endif

class gzip_stream {
    std::ofstream *file_stream;
    boost::iostreams::filtering_ostream *file_gzip_stream;
    
  public:
    gzip_stream(std::string file) {
          file_stream = new std::ofstream(file.c_str(), std::ios_base::binary);
          file_gzip_stream = new boost::iostreams::filtering_ostream;
          file_gzip_stream->push( boost::iostreams::gzip_compressor() );
          file_gzip_stream->push( *file_stream );
    }
    ~gzip_stream() {
        //delete file_gzip_stream;  
        //delete file_stream;
    }
    
    boost::iostreams::filtering_ostream *stream() {
        return file_gzip_stream;
    }
};

struct GlobalParams {

#ifdef USEMAGMA
    GPU_data gpu;
#endif
    
    int chain_num;
    
    int nRegions;
    int nLoci;
    int nInd;
    
    Rcpp::IntegerVector         nAlleles;
    
    arma::imat                  genotypes;
    Rcpp::IntegerVector         locate_indivs;
    arma::imat                  locate_genotypes;
    
    Rcpp::IntegerVector         indRegion;
    
    arma::mat                   pred_dist;
    
    std::vector<arma::mat>      count;          // NA  RxA[l]
    std::vector<arma::colvec>   sumCount;       // NA  Rx1
    std::vector<arma::mat>      allele_freq;    // NA  RxA[l]
    
    std::vector<arma::rowvec>   eta;            // NA  1xA[l]
    std::vector<double>         beta;           // 1
    std::vector<double>         xi;             // 1
    
    std::vector<double>         alpha;          // 4
    
    std::vector<arma::mat>      theta;          // NA  RxA[l]
    
    double                      anisoRatio;
    double                      anisoAngle;
    
    arma::mat                   anisoRatioMat;
    arma::mat                   anisoAngleMat;
    
    
    arma::mat                   dist;           // RxR
    arma::mat                   L;              // RxR
    arma::mat                   S;              // RxR
    arma::mat                   Sinv;           // RxR
    double                      Sdet;
    
    std::vector<double>         logLik_theta;   // NA  RxA[l]
    std::vector<double>         logLik_f;       // NA  RxA[l]
    
    

    arma::mat                   locs;           // Rx2
    arma::mat                   locsTrans;      // Rx2

    
    arma::mat                   boundary;
    
    arma::rowvec                theta_sd;        // NA
    arma::rowvec                mu_sd;       
    arma::rowvec                eta_sd;      
    arma::rowvec                xi_sd;
    arma::rowvec                beta_sd;
    arma::rowvec                alpha_sd;
    double                      angle_sd;
    double                      ratio_sd;
    
    arma::urowvec                thetaAttempt;        // NA
    arma::urowvec                thetaAccept;         // NA
    
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
    
    //std::vector<gzip_stream> cvfileStreams;
    //std::vector<std::vector<gzip_stream> > alfileStreams;
    
    std::vector<boost::iostreams::filtering_ostream*> cvfileGzStreams;
    std::vector<std::ofstream*> cvfileStreams;
    std::vector<std::vector<boost::iostreams::filtering_ostream*> > alfileGzStreams;
    std::vector<std::vector<std::ofstream*> > alfileStreams;
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
    
    bool OUTPUTALFREQ;
    bool GZIPOUTPUT;
        
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

    double THETASD;
    
    double LOCALSD;
    double GLOBALSD;
    
    double NULLPROB;
    double DELTA;
};

#endif