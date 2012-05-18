#include <Rmath.h>
#include <boost/lexical_cast.hpp>
#include <RcppArmadillo.h>

#include "Rscat_init.h"
#include "Rscat_cov.h"
#include "Rscat_util.h"

void calc_counts(GlobalParams &p, GlobalOptions &opt) {
    
    p.count.resize(p.nLoci);          // NA  RxA[l]
    p.sumCount.resize(p.nLoci);       // NA  Rx1 
    p.allele_freq.resize(p.nLoci);    // NA  RxA[l]
    
    for(int l=0; l<p.nLoci; ++l) {
        p.count[l] = arma::zeros<arma::mat>(p.nRegions,p.nAlleles[l]);
        
        if (opt.PSEUDOCOUNT != 0) 
            p.count[l].fill(opt.PSEUDOCOUNT);
        
        for(int i=0; i<p.nInd; i++) {
            int region = p.indRegion[i];
            
            for(int j=0; j<2; j++) {
                int allele = p.genotypes(2*i+j,l);
                
                if (allele >= 0)
                    p.count[l](region,allele)++;
            }
        }
        p.sumCount[l] = arma::sum(p.count[l],1);
        
        p.allele_freq[l] = p.count[l] / (p.sumCount[l] * arma::ones<arma::rowvec>(p.nAlleles[l]));
    }
}



void init_params(GlobalParams &p, GlobalOptions &opt) {
    
    p.alpha.resize(4);
    p.xi.resize(p.nLoci);
    p.eta.resize(p.nLoci);
    p.beta.resize(p.nLoci);
    
    p.theta.resize(p.nLoci);       
    p.logLik_theta.resize(p.nLoci);
    p.logLik_f.resize(p.nLoci);
        
    //p.alpha = std::vector<double>(4);
    //
    //p.eta  = std::vector<arma::rowvec>(p.nLoci);
    //p.beta = std::vector<double>(p.nLoci);
    //p.xi   = std::vector<double>(p.nLoci);
    //
    //p.theta = std::vector<arma::mat>(p.nLoci);
    //
    //p.logLik_theta = std::vector<double>(p.nLoci);
    //p.logLik_f     = std::vector<double>(p.nLoci);
    
    
    init_proposal_sd(p,opt);
    calc_counts(p,opt);
    
    p.anisoAngle = opt.ANGLE;
    p.anisoRatio = opt.RATIO;
    
    //for(int i=0; i<4; i++) {
    //    p.alpha[i] = (opt.FIXALPHA[i]) ? opt.ALPHA[i] : runif(1,opt.ALPHAMIN[i],opt.ALPHAMAX[i])[0]; 
    //}
    
    
    p.alpha[0] = (opt.FIXALPHA[0]) ? opt.ALPHA[0] : Rcpp::rgamma(1,2,1)[0];
    p.alpha[1] = (opt.FIXALPHA[1]) ? opt.ALPHA[1] : Rcpp::runif(1,opt.ALPHAMIN[1],500)[0]; 
    p.alpha[2] = (opt.FIXALPHA[2]) ? opt.ALPHA[2] : Rcpp::runif(1,opt.ALPHAMIN[2],opt.ALPHAMAX[2])[0]; 
    p.alpha[3] = (opt.FIXALPHA[3]) ? opt.ALPHA[3] : Rcpp::rgamma(1,2,0.2)[0];
    
    
    if (opt.FIXXI) {
        if (opt.XI.size() == 1)
            opt.XI = Rcpp::NumericVector(p.nLoci, opt.XI[0]);
        
        if (opt.XI.size() != p.nLoci)
            Rf_error("Length of values for fixed xi does not match number of Loci.");
    }
    if (opt.FIXBETA) {
        if (opt.BETA.size() == 1)
            opt.BETA = Rcpp::NumericVector(p.nLoci, opt.BETA[0]);
        
        if (opt.BETA.size() != p.nLoci)
            Rf_error("Length of values for fixed beta does not match number of Loci.");
    }
    if (opt.FIXETA) {
        if (opt.ETA.size() == 1)
            opt.ETA = Rcpp::NumericVector(p.nLoci, opt.ETA[0]);
        
        if (opt.ETA.size() != p.nLoci)
            Rf_error("Length of values for fixed eta does not match number of Loci.");
    }
    
    for(int l=0; l<p.nLoci; ++l) {
        
        p.xi[l] = (opt.FIXXI) ? opt.XI[l] 
                              : Rcpp::rnorm(1,0,0.25)[0]; //runif(1,opt.XIRANGE[0],opt.XIRANGE[1])[0];
        
        p.beta[l] = (opt.FIXBETA) ? opt.BETA[l]
                                  : abs(Rcpp::rnorm(1,3,1)[0]); //runif(1,opt.BETARANGE[0],opt.BETARANGE[1])[0]; 
        if (opt.FIXETA) {
            p.eta[l] = opt.ETA[l] * arma::ones<arma::rowvec>(p.nAlleles[l]);
        } else {
            p.eta[l] = p.beta[l]  * arma::randn<arma::rowvec>(p.nAlleles[l]);
        }
    }
    
    init_attempts(p);
    calc_params(p,opt);
}

void calc_params(GlobalParams &p, GlobalOptions &opt) {
    
    p.anisoAngleMat = calc_rotation_mat(p.anisoAngle);
    p.anisoRatioMat = calc_stretch_mat( p.anisoRatio);
    
    p.locsTrans = p.locs * p.anisoAngleMat * p.anisoRatioMat;
    
    p.dist = distance_mat(p.locsTrans);
    p.S    = calc_Sigma(p.alpha, p.dist, opt.USEMATERN);
    p.L    = calc_L(p.S);
    p.Sdet = det(p.S);
    
    arma::mat Linv = arma::inv( arma::trimatl(p.L) );
    p.Sinv = Linv.t() * Linv;

    for(int l=0; l<p.nLoci; ++l) {
        arma::mat mean = arma::ones<arma::colvec>(p.nRegions) * (p.xi[l] * p.eta[l]);
        p.theta[l] = mean +  p.L * arma::randn<arma::mat>(p.nRegions,p.nAlleles[l]);
        
        p.logLik_theta[l] = calc_multivar_normal_loglik(p.theta[l], mean, p.Sinv, p.Sdet);
        p.logLik_f[l] = calc_multinomial_loglik(p.theta[l], p.count[l], p.sumCount[l]);
    }
}

void init_proposal_sd(GlobalParams &p, GlobalOptions &opt) {
    
    p.theta_sd = arma::rowvec(p.nLoci);
    p.theta_sd.fill(opt.THETASD);
    
    p.eta_sd = arma::rowvec(p.nLoci);
    p.eta_sd.fill(opt.ETASD);
    
    p.xi_sd = arma::rowvec(p.nLoci);
    p.xi_sd.fill(opt.XISD);    

    p.beta_sd = arma::rowvec(p.nLoci);
    p.beta_sd.fill(opt.BETASD);
    
    p.alpha_sd = Rcpp::as<arma::rowvec>(opt.ALPHASD);
    
    p.angle_sd = opt.ANGLESD;
    p.ratio_sd = opt.RATIOSD;
}

void init_attempts(GlobalParams &p) {
    
    p.thetaAttempt = arma::zeros<arma::urowvec>(p.nLoci);
    p.thetaAccept  = arma::zeros<arma::urowvec>(p.nLoci);
    
    p.etaAttempt = arma::zeros<arma::urowvec>(p.nLoci);
    p.etaAccept  = arma::zeros<arma::urowvec>(p.nLoci);
    
    p.betaAttempt = arma::zeros<arma::urowvec>(p.nLoci);
    p.betaAccept  = arma::zeros<arma::urowvec>(p.nLoci);
    
    p.xiAttempt = arma::zeros<arma::urowvec>(p.nLoci);
    p.xiAccept  = arma::zeros<arma::urowvec>(p.nLoci);
    
    p.alphaAttempt = arma::zeros<arma::urowvec>(p.alpha.size());
    p.alphaAccept  = arma::zeros<arma::urowvec>(p.alpha.size());
    
    p.ratioAttempt = 0;
    p.ratioAccept  = 0;
    
    p.angleAttempt = 0;
    p.angleAccept  = 0;
}



void parseOptions(SEXP sexpRopt, GlobalOptions &opt)
{
    Rcpp::List Ropt(sexpRopt);
    
    opt.VERBOSE      = Rcpp::as<bool>(Ropt["VERBOSE"]);
    
    opt.TMPDIR       = Rcpp::as<std::string>(Ropt["TMPDIR"]);
    opt.FILESUFFIX   = Rcpp::as<std::string>(Ropt["FILESUFFIX"]);
    
    opt.ADAPT        = Rcpp::as<bool>(Ropt["ADAPT"]);
    opt.TUNEINTERVAL = Rcpp::as<int>(Ropt["TUNEINTERVAL"]);
    
    opt.LOCATE       = Rcpp::as<bool>(Ropt["LOCATE"]);
    opt.MAXCELL      = Rcpp::as<double>(Ropt["MAXCELL"]);
    
    opt.OUTPUTALFREQ = Rcpp::as<bool>(Ropt["OUTPUTALFREQ"]);
    opt.GZIPOUTPUT   = Rcpp::as<bool>(Ropt["GZIPOUTPUT"]);
    
    opt.RETURNFIT    = Rcpp::as<bool>(Ropt["RETURNFIT"]);
    opt.USEMATERN    = Rcpp::as<bool>(Ropt["USEMATERN"]);
    
    opt.PSEUDOCOUNT  = Rcpp::as<double>(Ropt["PSEUDOCOUNT"]);
    
    opt.FIXALPHA     = Rcpp::as<Rcpp::LogicalVector>(Ropt["FIXALPHA"]);
    opt.ALPHA        = Rcpp::as<Rcpp::NumericVector>(Ropt["ALPHA"]);
    opt.ALPHAMIN     = Rcpp::as<Rcpp::NumericVector>(Ropt["ALPHAMIN"]);
    opt.ALPHAMAX     = Rcpp::as<Rcpp::NumericVector>(Ropt["ALPHAMAX"]);
    
    opt.ALPHASD      = Rcpp::as<Rcpp::NumericVector>(Ropt["ALPHASD"]);
    
    opt.ANGLE        = Rcpp::as<double>(Ropt["ANGLE"]);
    opt.FIXANGLE     = Rcpp::as<bool>(Ropt["FIXANGLE"]);
    opt.ANGLESD      = Rcpp::as<double>(Ropt["ANGLESD"]);
    
    opt.RATIO        = Rcpp::as<double>(Ropt["RATIO"]);
    opt.FIXRATIO     = Rcpp::as<bool>(Ropt["FIXRATIO"]);
    opt.RATIOSD      = Rcpp::as<double>(Ropt["RATIOSD"]);
    
    opt.XIRANGE      = Rcpp::as<Rcpp::NumericVector>(Ropt["XIRANGE"]);
    opt.FIXXI        = Rcpp::as<bool>(Ropt["FIXXI"]);
    opt.XI           = Rcpp::as<Rcpp::NumericVector>(Ropt["XI"]);
    opt.XISD         = Rcpp::as<double>(Ropt["XISD"]);
    opt.SIGMAXI      = Rcpp::as<double>(Ropt["SIGMAXI"]);
    
    opt.FIXETA       = Rcpp::as<bool>(Ropt["FIXETA"]);
    opt.ETA          = Rcpp::as<Rcpp::NumericVector>(Ropt["ETA"]);
    opt.ETASD        = Rcpp::as<double>(Ropt["ETASD"]);
    
    opt.BETARANGE    = Rcpp::as<Rcpp::NumericVector>(Ropt["BETARANGE"]);
    opt.FIXBETA      = Rcpp::as<bool>(Ropt["FIXBETA"]);
    opt.BETA         = Rcpp::as<Rcpp::NumericVector>(Ropt["BETA"]);
    opt.BETASD       = Rcpp::as<double>(Ropt["BETASD"]);
    opt.SIGMABETA    = Rcpp::as<double>(Ropt["SIGMABETA"]);
    
    
    opt.THETASD      = Rcpp::as<double>(Ropt["THETASD"]);
    
    opt.LOCALSD      = Rcpp::as<double>(Ropt["LOCALSD"]);
    opt.GLOBALSD     = Rcpp::as<double>(Ropt["GLOBALSD"]);
    
    opt.NULLPROB     = Rcpp::as<double>(Ropt["NULLPROB"]);
    opt.DELTA        = Rcpp::as<double>(Ropt["DELTA"]);
}