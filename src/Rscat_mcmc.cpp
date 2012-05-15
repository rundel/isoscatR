#include <iomanip>
#include <boost/lexical_cast.hpp> 

#include <RcppArmadillo.h>

#include "Rscat_progressbar.h"
#include "Rscat_mcmc.h"
#include "Rscat_cov.h"
#include "Rscat_util.h"
#include "Rscat_init.h"
#include "Rscat_locate.h"

double tuneScale_1d(double acc_rate) {

    double scale = 1.0;

    if (acc_rate<0.001) {
        scale = 0.1;
    } else if (acc_rate<0.05) {
        scale = 0.5;
    } else if (acc_rate<0.2) {
        scale = 0.9;
    } else if (acc_rate>0.5) {
        scale = 1.1;
    } else if (acc_rate>0.75) {
        scale = 2.0;
    } else if (acc_rate>0.95) {
        scale = 10.0;
    }  
    
    return(scale);
}

double tuneScale_md(double acc_rate) {

    double scale = 1.0;

    if (acc_rate<0.001) {
        scale = 0.1;
    } else if (acc_rate<0.05) {
        scale = 0.5;
    } else if (acc_rate<0.1) {
        scale = 0.9;
    } else if (acc_rate>0.3) {
        scale = 1.1;
    } else if (acc_rate>0.5) {
        scale = 2.0;
    } else if (acc_rate>0.75) {
        scale = 10.0;
    }  
    
    return(scale);
}




void MCMCTune(GlobalParams &p, GlobalOptions &opt) {

    for(int l=0; l<p.alpha.size(); l++) {
        p.alpha_sd(l) *= tuneScale_1d( accept_ratio(p.alphaAccept(l),p.alphaAttempt(l)) );
    }
    
    p.angle_sd *= tuneScale_1d( accept_ratio(p.angleAccept,p.angleAttempt) );
    p.ratio_sd *= tuneScale_1d( accept_ratio(p.ratioAccept,p.ratioAttempt) );
    
    for(int l=0; l<p.nLoci; ++l) {
        p.xi_sd(l)    *= tuneScale_1d( accept_ratio(p.xiAccept(l), p.xiAttempt(l)) );
        p.beta_sd(l)  *= tuneScale_1d( accept_ratio(p.betaAccept(l), p.betaAttempt(l)) );
        
        p.eta_sd(l)   *= tuneScale_md( accept_ratio(p.etaAccept(l), p.etaAttempt(l)) );
        
        p.theta_sd[l] *= tuneScale_md( accept_ratio(p.thetaAccept[l], p.thetaAttempt[l]) );
        
    }
    
    init_attempts(p);
}


Rcpp::List MCMCLoop(GlobalParams &p, GlobalOptions &opt, int Niter, int Nthin, bool burnin, bool tune) {

    arma::mat samp_alpha, samp_mu, samp_xi, samp_beta, samp_aniso;
    std::vector<arma::mat> samp_eta(p.nLoci);
    arma::mat theta_fit;

    std::string prefix = tune ? "Tuning" : (burnin ? "Burnin" : "Sampling");

    progress_display progress_bar(Niter, std::cout, prefix);

    if(!burnin) {
        
        samp_alpha = arma::mat(Niter,p.alpha.size());
        samp_xi = arma::mat(Niter,p.nLoci);
        samp_beta = arma::mat(Niter,p.nLoci);
        
        samp_eta.resize(p.nLoci);
        for(int l=0; l<p.nLoci; l++) {
            samp_eta[l] = arma::mat(Niter,p.nAlleles[l]);
        }
        samp_aniso = arma::mat(Niter,2);

        theta_fit = arma::mat(Niter,1 + p.nRegions);
    }

    for (int i=0; i< Niter; i++){

        for (int j=0; j < Nthin; j++){
            MCMCStep(p, opt, burnin);
        }
        
        if(burnin && tune && i*Nthin % opt.TUNEINTERVAL == 0) 
            MCMCTune(p,opt);
        
        if (!burnin) {
            
            if (opt.LOCATE)
                update_location(p, opt);
            
            samp_aniso(i,0) = p.anisoAngle;
            samp_aniso(i,1) = p.anisoRatio;
            
            samp_alpha(i,0) = p.alpha[0];
            samp_alpha(i,1) = p.alpha[1];
            samp_alpha(i,2) = p.alpha[2];
            samp_alpha(i,3) = p.alpha[3];

            for (int l=0; l < p.nLoci; l++) {
                samp_xi(i,l)   = p.xi[l];
                samp_beta(i,l) = p.beta[l];
                
                samp_eta[l].row(i) = p.eta[l];
                
                //deviance.row(i) += -2*trans(p.logLik[l]);
            }
            
            if (opt.RETURNFIT) {
                double total_ss = 0;
                arma::colvec region_ss = arma::zeros<arma::colvec>(p.nRegions); 
                
                for(int l=0; l < p.nLoci; l++) {
                    arma::colvec rss = sum(square(p.allele_freq[l] - calc_f(p.theta[l])),1);
                    region_ss += rss;
                    total_ss += arma::accu(rss); 
                }
                
                theta_fit(i,0) = total_ss;
                for(int r=0; r<p.nRegions; ++r) 
                    theta_fit(i,r+1) = region_ss[r];
            }
        }
        
        ++progress_bar;
    }
    
    Rcpp::List res;
    
    if (!burnin) {
        
        std::vector<std::string> names_alpha, names_aniso, names_xibeta;
        names_alpha.push_back("alpha[0] - sigma2");
        names_alpha.push_back("alpha[1] - phi");
        names_alpha.push_back("alpha[2] - nu");
        names_alpha.push_back("alpha[3] - tau");
        
        names_aniso.push_back("anisotropy angle");
        names_aniso.push_back("anisotropy ratio");
        
        for (int l=0; l < p.nLoci; l++)
            names_xibeta.push_back("|xi|*beta [" + boost::lexical_cast<std::string>(l) + "]");
        
        
        res["alpha"] = Rcpp::List::create(Rcpp::Named("names") = names_alpha,
                                          Rcpp::Named("values") = samp_alpha);
                                    
        res["xibeta"] = Rcpp::List::create(Rcpp::Named("names") = names_xibeta,
                                           Rcpp::Named("values") = abs(samp_xi) % samp_beta);
        
        //res["deviance"] = List::create(Named("names") = "deviance",
        //                               Named("values") = sum(deviance,1));
        
        if (!opt.FIXRATIO && !opt.FIXANGLE) {
            res["aniso"] = Rcpp::List::create(Rcpp::Named("names") = names_aniso,
                                              Rcpp::Named("values") = samp_aniso);
        }
                                 
        if (opt.RETURNFIT) {
            
            std::vector<std::string> names_theta;
            names_theta.push_back("theta rss - total");
            
            for(int r=0; r<p.nRegions; ++r)
                names_theta.push_back("theta rss - Reg " + boost::lexical_cast<std::string>(r) ); 
            
            
            res["theta"] = Rcpp::List::create(Rcpp::Named("names") = names_theta,
                                              Rcpp::Named("values") = theta_fit);
        }
        
        // Calculating DIC based on fomula 6.12 on pg 183 of Bayesian Data Analysis
        //
        //std::vector<double> mean_alpha(p.alpha.size());
        //for(int i=0; i<p.alpha.size(); i++) {
        //    mean_alpha[i] = mean(samp_alpha.col(i));
        //}
        //
        //arma::mat mean_S = calc_Sigma(mean_alpha, p.dist, opt.USEMATERN);
        //arma::mat mean_L = calc_L(mean_S);
        //
        //arma::colvec D_theta_hat = zeros<arma::colvec>(p.nRegions);
        //for (int l=0; l < p.nLoci; l++) {
        //    double mean_xi = mean(samp_xi.col(l));
        //    
        //    arma::rowvec mean_eta = mean(samp_eta[l],0);
        //    arma::mat mean_X(p.nRegions,p.nAlleles[l]);
        //    
        //    for(int a=0; a<p.nAlleles[l]; a++) {
        //        mean_X.col(a) = trans(mean(samp_X[l][a],0));
        //    }        
        //    arma::mat mean_theta = arma::ones<arma::colvec>(p.nRegions) * (mean_xi * mean_eta) + mean_L * mean_X;
        //
        //    D_theta_hat += -2*calc_multinomial_loglik(mean_theta, p.count[l], p.sumCount[l]);
        //}
        //
        //
        //arma::colvec D_hat = trans(mean(deviance,0));
        //arma::colvec pV = trans(var(deviance,0,0)/2);
        //arma::colvec pD = D_hat - D_theta_hat;
        //
        //std::cout << "D_hat by region: ";
        //std::cout << trans(D_hat);
        //std::cout << "D_theta_hat by region: ";
        //std::cout << trans(D_theta_hat);
        
        //if (opt.VERBOSE) {
        //    std::cout << "Deviance Results:" << endl;
        //    std::cout << "=============================================" << endl;
        //
        //    std::cout << "D hat: " << setprecision(6) << setw(8) << arma::accu(D_hat) << endl;
        //    std::cout << "pD: " << setprecision(6) << setw(8) << arma::accu(pD) << endl;
        //    std::cout << "pV: " << setprecision(6) << setw(8) << arma::accu(pV) << endl;
        //
        //    std::cout << "DIC total (pD): " << setprecision(6) << setw(8) << arma::accu(D_hat+pD) << endl;
        //    std::cout << "DIC by region:";
        //    for(int r=0; r<p.nRegions; r++)
        //        std::cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pD(r)+0.5);
        //    std::cout << endl;
        //
        //    std::cout << "DIC total (pV): " << setprecision(6) << setw(8) << arma::accu(D_hat+pV) << endl;
        //    std::cout << "DIC by region:";
        //    for(int r=0; r<p.nRegions; r++)
        //        std::cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pV(r)+0.5);
        //    std::cout << endl << endl;
        //}
        
    }
    return(res);
}


void MCMCStep( GlobalParams &p, GlobalOptions &opt, bool burnin) {
    
    update_beta(p, opt);    //std::cout << "beta" << endl;
    update_eta(p,opt);      //std::cout << "eta" << endl;
    update_xi(p,opt);       //std::cout << "xi" << endl;
    update_alpha(p, opt);   //std::cout << "alpha" << endl;
    update_theta(p, opt);       //std::cout << "X" << endl;

    //update_anisotropy(p,opt);
}

/*
void update_anisotropy(GlobalParams &p, GlobalOptions &opt) {
    
    if (!opt.FIXRATIO) {
        p.ratioAttempt += 1;
        
        double newRatio;
        do {
            newRatio = p.anisoRatio + Rcpp::rnorm(1,0,p.ratio_sd)[0];
        } while (newRatio < 1);
        
        arma::mat newRatioMat = calc_stretch_arma::mat (newRatio);
        arma::mat newLocs = p.locs * p.anisoAngleMat * newRatioMat;
        arma::mat newDist = distance_arma::mat(newLocs);
        arma::mat newL    = calc_L(p.alpha, newDist, opt.USEMATERN);
        
        std::vector<arma::mat> newTheta(p.nLoci);
        std::vector<arma::colvec> newLogLik(p.nLoci);
        
        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.eta[l], p.xi[l], newL, p.X[l]);
            
            newLogLik[l] = calc_multinomial_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += arma::accu(newLogLik[l] - p.logLik[l]);
        }
        
        if( Rcpp::runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.ratioAccept +=1;
            
            p.anisoRatio   = newRatio;
            p.locsTrans    = newLocs;
            p.dist         = newDist;
            p.L            = newL;
            p.theta        = newTheta;
        }
    }
    
    if (!opt.FIXANGLE) {
        p.angleAttempt += 1;
        double newAngle = fmod(p.anisoAngle + Rcpp::rnorm(1,0,p.angle_sd)[0],  2*arma::math::pi());

        arma::mat newAngleMat = calc_rotation_arma::mat (newAngle);
        arma::mat newLocs = p.locs * newAngleMat * p.anisoRatioMat;
        arma::mat newDist = distance_arma::mat(newLocs);
        arma::mat newL    = calc_L(p.alpha, newDist, opt.USEMATERN);

        std::vector<arma::mat> newTheta(p.nLoci);
        std::vector<arma::colvec> newLogLik(p.nLoci);

        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.eta[l], p.xi[l], newL, p.X[l]);

            newLogLik[l] = calc_multinomial_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += arma::accu(newLogLik[l] - p.logLik[l]);
        }

        if( Rcpp::runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.angleAccept +=1;

            p.anisoAngle   = newAngle;
            p.locsTrans    = newLocs;
            p.dist         = newDist;
            p.L            = newL;
            p.theta        = newTheta;
        }
    }
}
*/


void update_alpha(GlobalParams &p, GlobalOptions &opt) {   // 4              
    

    for(int i = 0; i < 4; i++){
        if (opt.FIXALPHA[i]) continue;
        
        p.alphaAttempt[i] +=1;
        
        std::vector<double> newAlpha(p.alpha);
        arma::mat newL, newSigma;
        do {
            double h = Rcpp::rnorm(1,0,p.alpha_sd(i))[0];
            
            if (i == 1)
                newAlpha[i] = p.alpha[i] * exp( h );
            else
                newAlpha[i] = p.alpha[i] + h;
            
            
            if (newAlpha[i] < opt.ALPHAMIN[i] || newAlpha[i] > opt.ALPHAMAX[i])
                continue;
            
            newSigma = calc_Sigma(newAlpha, p.dist, opt.USEMATERN);
            newL = calc_L(newSigma);
            
        } while (newL.is_empty());
        
        arma::mat newLInv = arma::inv( arma::trimatl(newL) );
        arma::mat newSigmaInv = newLInv.t() * newLInv;
        double newSigmaDet = arma::det(newSigma);
        
        std::vector<double> newLogLik_theta(p.nLoci);
        
        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            
            arma::mat mean = arma::ones<arma::colvec>(p.nRegions) * (p.xi[l] * p.eta[l]);
            newLogLik_theta[l] = calc_multivar_normal_loglik(p.theta[l], mean, newSigmaInv, newSigmaDet);            
            
            logLikRatio += newLogLik_theta[l]-p.logLik_theta[l];
        }
        
        // FIXME - transition probabilities are not symmetric 
        
        if( Rcpp::runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.alphaAccept[i] +=1;
            
            p.alpha[i]     = newAlpha[i];
            p.L            = newL;
            p.Sinv         = newSigmaInv;
            p.Sdet         = newSigmaDet;
            p.logLik_theta = newLogLik_theta;
        }
        
    }
}

void update_theta(GlobalParams &p, GlobalOptions &opt) {
    
    // Based on Eqn A.2 from BDA
    
    for(int l=0; l<p.nLoci; l++) {
        
        p.thetaAttempt[l] += 1;
        
        arma::mat mean = arma::ones<arma::colvec>(p.nRegions) * (p.xi[l] * p.eta[l]);
        arma::mat norm = arma::randn<arma::mat>(p.nRegions,p.nAlleles[l]);
        
        arma::mat newTheta = p.theta[l] + p.theta_sd[l]*norm;
        
        double newLogLik_theta = calc_multivar_normal_loglik(newTheta, mean, p.Sinv, p.Sdet);
        double newLogLik_f     = calc_multinomial_loglik(newTheta, p.count[l], p.sumCount[l]);
        
        double logLikRatio = newLogLik_theta-p.logLik_theta[l]+newLogLik_f-p.logLik_f[l];
        
        if( Rcpp::runif(1)[0] < exp( logLikRatio ) ){
            p.thetaAccept[l]  +=1;
        
            p.theta[l]        = newTheta;
            p.logLik_theta[l] = newLogLik_theta;
            p.logLik_f[l]     = newLogLik_f;
        }
    }
}


void update_xi(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXXI)
        return;
    
    for(int l=0; l < p.nLoci; l++) {
        
        p.xiAttempt[l] +=1;
        
        double newXi = p.xi[l] + Rcpp::rnorm(1,0,p.xi_sd(l))[0];
        
        arma::mat new_mean = arma::ones<arma::colvec>(p.nRegions) * (newXi   * p.eta[l]);
        
        double newLogLik_theta = calc_multivar_normal_loglik(p.theta[l], new_mean, p.Sinv, p.Sdet);
        
        double logLikRatio = newLogLik_theta-p.logLik_theta[l];
        
        if( Rcpp::runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.xiAccept[l] += 1;

            p.xi[l]           = newXi;                
            p.logLik_theta[l] = newLogLik_theta;
        }
    }
}



void update_eta(GlobalParams &p, GlobalOptions &opt) {    // NA  1xA[l]
    
    if (opt.FIXETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        p.etaAttempt[l] += 1;
        
        arma::rowvec newEta = p.eta[l] + p.eta_sd(l) * arma::randn<arma::rowvec>(p.nAlleles[l]);
        
        arma::mat new_mean = arma::ones<arma::colvec>(p.nRegions) * (p.xi[l] * newEta  );
        
        double newLogLik_theta = calc_multivar_normal_loglik(p.theta[l], new_mean, p.Sinv, p.Sdet);
        
        //prior on eta is N(0,beta) 
        double logLikRatio = newLogLik_theta - p.logLik_theta[l] - 0.5*arma::accu(square(newEta) - square(p.eta[l])) / (p.beta[l]*p.beta[l]);
        
        if( Rcpp::runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.etaAccept[l] += 1;

            p.logLik_theta[l] = newLogLik_theta;
            p.eta[l] = newEta;                
        }
    }
}

// beta is the prior precision for Mu (ie Mu is N(0,Beta)
// prior on beta is Uniform(BETAMIN,BETAMAX)

// Latex Model Eqns:
// &\left( \prod_l \prod_{a[l]} \frac{1}{\sqrt{2\pi}\beta}e^{-\frac{1}{2}\frac{\eta_{la}^2}{\beta^2}} \right) p(\beta)\\
// &\left( \sum_l \sum_{a[l]} -\frac{1}{2}\log (2\pi) - \log \beta - \frac{1}{2} \frac{\eta_{la}^2}{\beta^2} \right)+\log(p(\beta))\\
// &-\left(\sum_l n_{a[l]}\right) \log \beta - \frac{1}{2\beta^2} \sum_l \sum_{a[l]}\eta_{la}^2+\log(p(\beta))

void update_beta(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXBETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        
        p.betaAttempt[l] += 1;
        
        double newBeta;
        do {
            newBeta = p.beta[l] + Rcpp::rnorm(1,0,p.beta_sd(l))[0];
        } while(newBeta < opt.BETARANGE[0] || newBeta > opt.BETARANGE[1]);
        
        double etasq_sum = arma::accu(p.eta[l] % p.eta[l]);
        
        double logLikRatio = -p.nAlleles[l]*(log(newBeta/p.beta[l])) - 0.5 * etasq_sum * (1/(newBeta*newBeta) - 1/(p.beta[l]*p.beta[l]));
        
        if( Rcpp::runif(1)[0] < exp( logLikRatio ) ){ //accept move 
            p.betaAccept[l] += 1;
            
            p.beta[l] = newBeta;
        }
    }
}
