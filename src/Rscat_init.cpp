#include <Rmath.h>
#include <boost/lexical_cast.hpp>
#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

void calc_counts(GlobalParams &p, GlobalOptions &opt) {
    
    p.count.resize(p.nLoci);          // NA  RxA[l]
    p.sumCount.resize(p.nLoci);       // NA  Rx1 
    p.allele_freq.resize(p.nLoci);    // NA  RxA[l]
    
    for(int l=0; l<p.nLoci; ++l) {
        p.count[l] = zeros<mat>(p.nRegions,p.nAlleles[l]);
        
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
        p.sumCount[l] = sum(p.count[l],1);
        
        p.allele_freq[l] = p.count[l] / (p.sumCount[l] * ones<rowvec>(p.nAlleles[l]));
    }
}



void init_params(GlobalParams &p, GlobalOptions &opt) {
    
    p.alpha.resize(4);
    p.xi.resize(p.nLoci);
    p.eta.resize(p.nLoci);
    p.beta.resize(p.nLoci);
    
    p.X.resize(p.nLoci);              // NA  RxA[l]
    p.theta.resize(p.nLoci);          // NA  RxA[l]
    p.logLik_theta.resize(p.nLoci);         // NA  Rx1
    p.logLik_f.resize(p.nLoci);         // NA  Rx1
    
    init_proposal_sd(p,opt);
    calc_counts(p,opt);
    
    p.anisoAngle = opt.ANGLE;
    p.anisoRatio = opt.RATIO;
    
    //for(int i=0; i<4; i++) {
    //    p.alpha[i] = (opt.FIXALPHA[i]) ? opt.ALPHA[i] : runif(1,opt.ALPHAMIN[i],opt.ALPHAMAX[i])[0]; 
    //}
    
    
    p.alpha[0] = (opt.FIXALPHA[0]) ? opt.ALPHA[0] : rgamma(1,2,1)[0];
    p.alpha[1] = (opt.FIXALPHA[1]) ? opt.ALPHA[1] : runif(1,opt.ALPHAMIN[1],500)[0]; 
    p.alpha[2] = (opt.FIXALPHA[2]) ? opt.ALPHA[2] : runif(1,opt.ALPHAMIN[2],opt.ALPHAMAX[2])[0]; 
    p.alpha[3] = (opt.FIXALPHA[3]) ? opt.ALPHA[3] : rgamma(1,2,0.2)[0];
    
    
    if (opt.FIXXI) {
        if (opt.XI.size() == 1)
            opt.XI = NumericVector(p.nLoci, opt.XI[0]);
        
        if (opt.XI.size() != p.nLoci)
            Rf_error("Length of values for fixed xi does not match number of Loci.");
    }
    if (opt.FIXBETA) {
        if (opt.BETA.size() == 1)
            opt.BETA = NumericVector(p.nLoci, opt.BETA[0]);
        
        if (opt.BETA.size() != p.nLoci)
            Rf_error("Length of values for fixed beta does not match number of Loci.");
    }
    if (opt.FIXETA) {
        if (opt.ETA.size() == 1)
            opt.ETA = NumericVector(p.nLoci, opt.ETA[0]);
        
        if (opt.ETA.size() != p.nLoci)
            Rf_error("Length of values for fixed eta does not match number of Loci.");
    }
    
    for(int l=0; l<p.nLoci; ++l) {
        
        p.xi[l] = (opt.FIXXI) ? opt.XI[l] 
                                : rnorm(1,0,0.25)[0]; //runif(1,opt.XIRANGE[0],opt.XIRANGE[1])[0];
        
        p.beta[l] = (opt.FIXBETA) ? opt.BETA[l]
                                    : abs(rnorm(1,3,1)[0]); //runif(1,opt.BETARANGE[0],opt.BETARANGE[1])[0]; 
        
        p.eta[l] = (opt.FIXETA) ? p.eta[l].fill(opt.ETA[l])
                                  : p.beta[l] * randn<rowvec>(p.nAlleles[l]);
        
        p.X[l] =  randn<mat>(p.nRegions,p.nAlleles[l]);
    }
    
    init_attempts(p);
    calc_params(p,opt);
}

void calc_params(GlobalParams &p, GlobalOptions &opt) {
    
    p.anisoAngleMat = calc_rotation_mat(p.anisoAngle);
    p.anisoRatioMat = calc_stretch_mat( p.anisoRatio);
    
    p.locsTrans = p.locs * p.anisoAngleMat * p.anisoRatioMat;
    
    p.dist = calc_distance_mat(p.locsTrans);
    p.S    = calc_Sigma(p.alpha, p.dist, opt.USEMATERN);
    p.L    = calc_L(p.S);
    p.Sdet = det(p.S);
    
    mat Linv = inv( trimatl(p.L) );
    p.Sinv = Linv.t() * Linv;

    for(int l=0; l<p.nLoci; ++l) {
        mat mean = ones<colvec>(p.nRegions) * (p.xi[l] * p.eta[l]);
        p.theta[l] = mean +  p.L * p.X[l];
        
        p.logLik_theta[l] = calc_multivar_normal_loglik(p.theta[l], mean, p.Sinv, p.Sdet);
        p.logLik_f[l] = calc_multinomial_loglik(p.theta[l], p.count[l], p.sumCount[l]);
    }
}

void init_proposal_sd(GlobalParams &p, GlobalOptions &opt) {
    
    
    p.x_sd.resize(p.nLoci);
    for(int l=0; l<p.nLoci; ++l) {
        p.x_sd[l] = colvec(p.nRegions);
        p.x_sd[l].fill(opt.XSD);        
    }
    
    p.eta_sd = rowvec(p.nLoci);
    p.eta_sd.fill(opt.ETASD);
    
    p.xi_sd = rowvec(p.nLoci);
    p.xi_sd.fill(opt.XISD);    

    p.beta_sd = rowvec(p.nLoci);
    p.beta_sd.fill(opt.BETASD);
    
    p.alpha_sd = as<rowvec>(opt.ALPHASD);
    
    p.angle_sd = opt.ANGLESD;
    p.ratio_sd = opt.RATIOSD;
}

void init_attempts(GlobalParams &p) {
    
    p.Xattempt.resize(p.nLoci);
    p.Xaccept.resize(p.nLoci);
    
    for(int l=0; l<p.nLoci; ++l) {
        p.Xattempt[l]  = zeros<ucolvec>(p.nRegions);
        p.Xaccept[l]   = zeros<ucolvec>(p.nRegions);
    }
    
    p.etaAttempt = zeros<urowvec>(p.nLoci);
    p.etaAccept  = zeros<urowvec>(p.nLoci);
    
    p.betaAttempt = zeros<urowvec>(p.nLoci);
    p.betaAccept  = zeros<urowvec>(p.nLoci);
    
    p.xiAttempt = zeros<urowvec>(p.nLoci);
    p.xiAccept  = zeros<urowvec>(p.nLoci);
    
    p.alphaAttempt = zeros<urowvec>(p.alpha.size());
    p.alphaAccept  = zeros<urowvec>(p.alpha.size());
    
    p.ratioAttempt = 0;
    p.ratioAccept  = 0;
    
    p.angleAttempt = 0;
    p.angleAccept  = 0;
}



void open_cvfiles(GlobalParams &p, GlobalOptions &opt) {
    
    p.cvfileStreams.resize(p.locate_indivs.size());
    p.cvfileGzStreams.resize(p.locate_indivs.size());

    for(int l=0; l<p.locate_indivs.size(); l++) {
        
        string file = opt.TMPDIR + "/CV_Ind" 
                      + boost::lexical_cast<string>(p.locate_indivs[l]) + "_"
                      + boost::lexical_cast<string>(p.chain_num) + ".gz";
        
        p.cvfileStreams[l] = new ofstream(file.c_str(), ios_base::binary);
        
        p.cvfileGzStreams[l] = new boost::iostreams::filtering_ostream;
        p.cvfileGzStreams[l]->push( boost::iostreams::gzip_compressor() );
        p.cvfileGzStreams[l]->push( *(p.cvfileStreams[l]) );
    
    }
}

void close_cvfiles(GlobalParams &p, GlobalOptions &opt) {
    
    for(int l=0; l<p.locate_indivs.size(); l++) {
        delete p.cvfileGzStreams[l];  
        delete p.cvfileStreams[l];
    }
}


void open_allelefiles(GlobalParams &p, GlobalOptions &opt) {
    
    p.alfileStreams.resize(p.nLoci);
    p.alfileGzStreams.resize(p.nLoci);
 
    for(int l=0; l<p.nLoci; l++) {
        
        p.alfileStreams[l].resize(p.nAlleles[l]);
        p.alfileGzStreams[l].resize(p.nAlleles[l]);
        
        for(int j=0; j<p.nAlleles[l]; j++) {
            string file = opt.TMPDIR + "/Al" 
                          + boost::lexical_cast<string>(l+1) + "-"
                          + boost::lexical_cast<string>(j+1) + "_"
                          + boost::lexical_cast<string>(p.chain_num) + ".gz";
        
            p.alfileStreams[l][j] = new ofstream(file.c_str(), ios_base::binary);
        
            p.alfileGzStreams[l][j] = new boost::iostreams::filtering_ostream;
            p.alfileGzStreams[l][j]->push( boost::iostreams::gzip_compressor() );
            p.alfileGzStreams[l][j]->push( *(p.alfileStreams[l][j]) );
        }
    }
}
 
void close_allelefiles(GlobalParams &p, GlobalOptions &opt) {
    
    for(int l=0; l<p.nLoci; l++) {
        for(int j=0; j<p.nAlleles[l]; j++) {
            delete p.alfileGzStreams[l][j];  
            delete p.alfileStreams[l][j];
        }
    }
}