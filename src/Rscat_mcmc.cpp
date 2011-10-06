#include <vector>
#include <sstream>
#include "Rscat.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

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
        p.alpha_sd(l) *= tuneScale_1d( calc_accept_ratio(p.alphaAccept(l),p.alphaAttempt(l)) );
    }
    
    p.angle_sd *= tuneScale_1d( calc_accept_ratio(p.angleAccept,p.angleAttempt) );
    p.ratio_sd *= tuneScale_1d( calc_accept_ratio(p.ratioAccept,p.ratioAttempt) );
    
    for(int l=0; l<p.nLoci; ++l) {
        p.mu_sd(l)    *= tuneScale_1d( calc_accept_ratio(p.muAccept(l), p.muAttempt(l)) );
        p.xi_sd(l)    *= tuneScale_1d( calc_accept_ratio(p.xiAccept(l), p.xiAttempt(l)) );
        p.beta_sd(l)  *= tuneScale_1d( calc_accept_ratio(p.betaAccept(l), p.betaAttempt(l)) );
        
        p.eta_sd(l)   *= tuneScale_md( calc_accept_ratio(p.etaAccept(l), p.etaAttempt(l)) );
        
        for(int r=0; r<p.nRegions; ++r) {
            p.x_sd[l](r) *= tuneScale_md( calc_accept_ratio(p.Xaccept[l](r), p.Xattempt[l](r)) );
        }
    }
    
    init_attempts(p);

}


List MCMCLoop(GlobalParams &p, GlobalOptions &opt, int Niter, int Nthin, bool burnin, bool tune) {

    mat samp_alpha, samp_mu, samp_xi, samp_beta;
    vector<mat> samp_eta(p.nLoci);
    vector< vector<mat> > samp_X(p.nLoci);
    mat samp_aniso, deviance, theta_fit;

    if(!burnin) {
        
        samp_alpha = mat(Niter,p.alpha.size());
        samp_mu = mat(Niter,p.nLoci);
        samp_xi = mat(Niter,p.nLoci);
        samp_beta = mat(Niter,p.nLoci);
        
        samp_eta.resize(p.nLoci);
        samp_X.resize(p.nLoci);
        for(int l=0; l<p.nLoci; l++) {
            samp_eta[l] = mat(Niter,p.nAlleles[l]);
            
            samp_X[l].resize(p.nAlleles[l]);
            for(int i=0; i<p.nAlleles[l]; i++) {
                samp_X[l][i] = mat(Niter,p.nRegions);
            }
        }
        theta_fit = mat(Niter,1 + p.nRegions);
        deviance = zeros<mat>(Niter,p.nRegions);
        samp_aniso = mat(Niter,2);
    }

    for (int i=0; i< Niter; i++){
        if (opt.VERBOSE && i*Nthin % 1000 == 0)
            cout << "Iter " << setw(6) << i*Nthin << endl;
        
        for (int j=0; j < Nthin; j++){
            MCMCStep(p, opt, burnin);
        }
        
        if(burnin && tune && i*Nthin % opt.TUNEINTERVAL == 0) 
            MCMCTune(p,opt);
        
        if (!burnin) {
            if (opt.LOCATE)
                update_Location(p,opt);
            if (opt.CROSSVALIDATE)
                update_LocationCV(p,opt);
            
            samp_aniso(i,0) = p.anisoAngle;
            samp_aniso(i,1) = p.anisoRatio;
            
            samp_alpha(i,0) = p.alpha[0];
            samp_alpha(i,1) = p.alpha[1];
            samp_alpha(i,2) = p.alpha[2];
            samp_alpha(i,3) = p.alpha[3];

            for (int l=0; l < p.nLoci; l++) {
                samp_mu(i,l)   = p.mu[l];
                samp_xi(i,l)   = p.xi[l];
                samp_beta(i,l) = p.beta[l];
                
                samp_eta[l].row(i) = p.eta[l];
                for(int a=0; a<p.nAlleles[l]; a++) {
                    samp_X[l][a].row(i) = trans(p.X[l].col(a));
                }
                
                deviance.row(i) += -2*trans(p.logLik[l]);
            }
            
            if (opt.RETURNFIT) {
                double total_ss = 0;
                colvec region_ss = zeros<colvec>(p.nRegions); 
                
                for(int l=0; l < p.nLoci; l++) {
                    colvec rss = sum(square(p.allele_freq[l] - calc_f(p.theta[l])),1);
                    region_ss += rss;
                    total_ss += accu(rss); 
                }
                
                theta_fit(i,0) = total_ss;
                for(int r=0; r<p.nRegions; ++r) 
                    theta_fit(i,r+1) = region_ss[r];
            }
        }
    }
    
    List res;
    
    if (!burnin) {
        //res = List::create(Named("alpha") = samp_alpha,
        //                   Named("xibeta") = samp_xibeta,
        //                   Named("mu") = samp_mu);
        
        vector<string> names_alpha, names_aniso, names_xibeta, names_mu;
        names_alpha.push_back("alpha[0] - sigma2");
        names_alpha.push_back("alpha[1] - phi");
        names_alpha.push_back("alpha[2] - nu");
        names_alpha.push_back("alpha[3] - tau");
        
        names_aniso.push_back("anisotropy angle");
        names_aniso.push_back("anisotropy ratio");
        
        for (int l=0; l < p.nLoci; l++) {
            stringstream ss1, ss2;
            
            ss1 << "|xi|*beta [" << l << "]";
            names_xibeta.push_back(ss1.str());
            
            ss2 << "mu["<< l << "]";
            names_mu.push_back(ss2.str());
        }
        
        
        res["alpha"] = List::create(Named("names") = names_alpha,
                                    Named("values") = samp_alpha);
                                    
        res["xibeta"] = List::create(Named("names") = names_xibeta,
                                     Named("values") = abs(samp_xi) % samp_beta);
        
        res["deviance"] = List::create(Named("names") = "deviance",
                                       Named("values") = sum(deviance,1));
        
        if (!opt.FIXMU) {
            res["mu"] = List::create(Named("names") = names_mu,
                                     Named("values") = samp_mu);
        }
        
        if (!opt.FIXRATIO && !opt.FIXANGLE) {
            res["aniso"] = List::create(Named("names") = names_aniso,
                                        Named("values") = samp_aniso);
        }
                                 
        if (opt.RETURNFIT) {
            
            vector<string> names_theta;
            names_theta.push_back("theta rss - total");
            
            for(int r=0; r<p.nRegions; ++r) {
                stringstream ss;
                ss << "theta rss - Reg " << r;
                names_theta.push_back(ss.str()); 
            }
            
            res["theta"] = List::create(Named("names") = names_theta,
                                        Named("values") = theta_fit);
        }
        
        // Calculating DIC based on fomula 6.12 on pg 183 of Bayesian Data Analysis
        
        vector<double> mean_alpha(p.alpha.size());
        for(int i=0; i<p.alpha.size(); i++) {
            mean_alpha[i] = mean(samp_alpha.col(i));
        }
        
        mat mean_L = calc_L(mean_alpha, p.dist, opt.USEMATERN);
        
        colvec D_theta_hat = zeros<colvec>(p.nRegions);
        for (int l=0; l < p.nLoci; l++) {
            double mean_mu = mean(samp_mu.col(l));
            double mean_xi = mean(samp_xi.col(l));
            
            rowvec mean_eta = mean(samp_eta[l],0);
            mat mean_X(p.nRegions,p.nAlleles[l]);
            
            for(int a=0; a<p.nAlleles[l]; a++) {
                mean_X.col(a) = trans(mean(samp_X[l][a],0));
            }        
            mat mean_theta = calc_theta(mean_mu, mean_eta, mean_xi, mean_L, mean_X);
    
            D_theta_hat += -2*calc_multinom_loglik(mean_theta, p.count[l], p.sumCount[l]);
        }
        
        
        colvec D_hat = trans(mean(deviance,0));
        colvec pV = trans(var(deviance,0,0)/2);
        colvec pD = D_hat - D_theta_hat;
        
        //cout << "D_hat by region: ";
        //cout << trans(D_hat);
        //cout << "D_theta_hat by region: ";
        //cout << trans(D_theta_hat);
        
        if (opt.VERBOSE) {
            cout << "Deviance Results:" << endl;
            cout << "=============================================" << endl;
        
            cout << "D hat: " << setprecision(6) << setw(8) << accu(D_hat) << endl;
            cout << "pD: " << setprecision(6) << setw(8) << accu(pD) << endl;
            cout << "pV: " << setprecision(6) << setw(8) << accu(pV) << endl;
        
            cout << "DIC total (pD): " << setprecision(6) << setw(8) << accu(D_hat+pD) << endl;
            cout << "DIC by region:";
            for(int r=0; r<p.nRegions; r++)
                cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pD(r)+0.5);
            cout << endl;
        
            cout << "DIC total (pV): " << setprecision(6) << setw(8) << accu(D_hat+pV) << endl;
            cout << "DIC by region:";
            for(int r=0; r<p.nRegions; r++)
                cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pV(r)+0.5);
            cout << endl << endl;
        }
        
    }
    return(res);
}


void MCMCStep( GlobalParams &p, GlobalOptions &opt, bool burnin) {
    
    update_Beta(p, opt);    //cout << "beta" << endl;
    update_Eta(p,opt);      //cout << "eta" << endl;
    update_Xi(p,opt);       //cout << "xi" << endl;
    update_Mu(p, opt);      //cout << "mu" << endl;
    update_Alpha(p, opt);   //cout << "alpha" << endl;
    update_X(p, opt);       //cout << "X" << endl;
    
    update_anisotropy(p,opt);
}

void update_anisotropy(GlobalParams &p, GlobalOptions &opt) {
    
    if (!opt.FIXRATIO) {
        p.ratioAttempt += 1;
        
        double newRatio;
        do {
            newRatio = p.anisoRatio + rnorm(1,0,p.ratio_sd)[0];
        } while (newRatio < 1);
        
        mat newRatioMat = calc_stretch_mat (newRatio);
        mat newLocs = p.locs * p.anisoAngleMat * newRatioMat;
        mat newDist = calc_distance_mat(newLocs);
        mat newL    = calc_L(p.alpha, newDist, opt.USEMATERN);
        
        vector<mat> newTheta(p.nLoci);
        vector<colvec> newLogLik(p.nLoci);
        
        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.mu[l], p.eta[l], p.xi[l], newL, p.X[l]);
            
            newLogLik[l] = calc_multinom_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += accu(newLogLik[l] - p.logLik[l]);
        }
        
        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.ratioAccept +=1;
            
            p.anisoRatio   = newRatio;
            p.locsTrans    = newLocs;
            p.dist         = newDist;
            p.L            = newL;
            p.theta        = newTheta;
            p.logLik       = newLogLik;
        }
    }
    
    if (!opt.FIXANGLE) {
        p.angleAttempt += 1;
        double newAngle = fmod(p.anisoAngle + rnorm(1,0,p.angle_sd)[0],  2*math::pi());

        mat newAngleMat = calc_rotation_mat (newAngle);
        mat newLocs = p.locs * newAngleMat * p.anisoRatioMat;
        mat newDist = calc_distance_mat(newLocs);
        mat newL    = calc_L(p.alpha, newDist, opt.USEMATERN);

        vector<mat> newTheta(p.nLoci);
        vector<colvec> newLogLik(p.nLoci);

        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.mu[l], p.eta[l], p.xi[l], newL, p.X[l]);

            newLogLik[l] = calc_multinom_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += accu(newLogLik[l] - p.logLik[l]);
        }

        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.angleAccept +=1;

            p.anisoAngle   = newAngle;
            p.locsTrans    = newLocs;
            p.dist         = newDist;
            p.L            = newL;
            p.theta        = newTheta;
            p.logLik       = newLogLik;
        }
    }
}



void update_Alpha(GlobalParams &p, GlobalOptions &opt) {   // 4              
    

    for(int i = 0; i < 4; i++){
        if (opt.FIXALPHA[i]) continue;
        
        p.alphaAttempt[i] +=1;
        
        vector<double> newAlpha(p.alpha);
        mat newL;
        do {
            double h = rnorm(1,0,p.alpha_sd(i))[0];
            
            switch(i) {
                case 1:
                    newAlpha[i] = p.alpha[i] * exp( h );
                    break;
                case 2:
                    newAlpha[i] = p.alpha[i] * exp( h );
                    break;
                default:
                    newAlpha[i] = p.alpha[i] + h;
            }
            
            
            if (newAlpha[i] < opt.ALPHAMIN[i] || newAlpha[i] > opt.ALPHAMAX[i])
                continue;
            
            newL = calc_L(newAlpha, p.dist,opt.USEMATERN);
            
        } while (newL.is_empty());
        
        // FIXME - transition probabilities are not symmetric 
        
        vector<mat> newTheta(p.nLoci);
        vector<colvec> newLogLik(p.nLoci);
        
        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.mu[l], p.eta[l], p.xi[l], newL, p.X[l]);
            
            newLogLik[l] = calc_multinom_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += accu(newLogLik[l] - p.logLik[l]);
        }
        
        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.alphaAccept[i] +=1;
            
            p.alpha[i]     = newAlpha[i];
            p.L            = newL;
            p.theta        = newTheta;
            p.logLik       = newLogLik;
        }
        
    }
}






// update Mu (the background "ancestral" allele freqs)
void update_Mu(GlobalParams &p, GlobalOptions &opt) {    // NA  1xA[l]
    
    if (opt.FIXMU) return;
    
    for(int l=0; l<p.nLoci; l++) {
        
        p.muAttempt[l] +=1;
        
        double newMu = p.mu[l] + rnorm(1,0,p.mu_sd(l))[0];
        
        mat newTheta = calc_theta(newMu, p.eta[l], p.xi[l], p.L, p.X[l]);
        
        colvec newLogLik = calc_multinom_loglik(newTheta,    p.count[l], p.sumCount[l]);
                                                           // Normal prior N(0,10)
        double logLikRatio = accu(newLogLik - p.logLik[l]) - 0.5*(pow(newMu,2)-pow(p.mu[l],2))/pow(10.0,2); 
        
        if( runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.muAccept[l] +=1;
            
            p.theta[l] = newTheta;
            p.logLik[l] = newLogLik;
            p.mu[l] = newMu;                
        }
    }
}

void update_Xi(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXXI) return;
    
    for(int l=0; l < p.nLoci; l++) {
        
        p.xiAttempt[l] +=1;
        
        double newXi = p.xi[l] + rnorm(1,0,p.xi_sd(l))[0];
        mat newTheta = calc_theta(p.mu[l], p.eta[l], newXi, p.L, p.X[l]);
        
        colvec newLogLik = calc_multinom_loglik(newTheta,    p.count[l], p.sumCount[l]);
        
        double logLikRatio = accu(newLogLik - p.logLik[l]); 
        
        if( runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.xiAccept[l] += 1;

            p.theta[l] = newTheta;
            p.logLik[l] = newLogLik;
            p.xi[l] = newXi;                
        }
    }
}




void update_Eta(GlobalParams &p, GlobalOptions &opt) {    // NA  1xA[l]
    
    if (opt.FIXETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        p.etaAttempt[l] += 1;
        
        rowvec newEta = p.eta[l] + p.eta_sd(l) * randn<rowvec>(p.nAlleles[l]);
        mat newTheta = calc_theta(p.mu[l], newEta, p.xi[l], p.L, p.X[l]);
        
        colvec newLogLik = calc_multinom_loglik(newTheta, p.count[l], p.sumCount[l]);
        
        //prior on eta is N(0,beta) 
        double logLikRatio = accu(newLogLik - p.logLik[l]) - 0.5*accu(square(newEta) - square(p.eta[l])) / (p.beta[l]*p.beta[l]);
        
        if( runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.etaAccept[l] += 1;

            p.theta[l] = newTheta;
            p.logLik[l] = newLogLik;
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

void update_Beta(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXBETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        
        p.betaAttempt[l] += 1;
        
        double newBeta;
        do {
            newBeta = p.beta[l] + rnorm(1,0,p.beta_sd(l))[0];
        } while(newBeta < opt.BETARANGE[0] || newBeta > opt.BETARANGE[1]);
        
        double etasq_sum = accu(square(p.eta[l]));
        
        double logLikRatio = -p.nAlleles[l]*(log(newBeta/p.beta[l])) - 0.5 * etasq_sum * (1/(newBeta*newBeta) - 1/(p.beta[l]*p.beta[l]));
        
        if( runif(1)[0] < exp( logLikRatio ) ){ //accept move 
            p.beta[l] = newBeta;
            p.betaAccept[l] += 1;
        }
    }
}



void update_X(GlobalParams &p, GlobalOptions &opt) {
    
    for(int l=0; l<p.nLoci; l++) {
        for(int r=0; r<p.nRegions; r++) {
            p.Xattempt[l](r) += 1;
        
            mat newX = p.X[l];
            newX.row(r) += p.x_sd[l](r) * randn<rowvec>(p.nAlleles[l]);
            
            //FIXME
            mat newTheta = calc_theta(p.mu[l], p.eta[l], p.xi[l], p.L, newX);

            colvec newLogLik = calc_multinom_loglik(newTheta,    p.count[l], p.sumCount[l]);
        
            double logLikRatio = accu(newLogLik - p.logLik[l]) - 0.5 * accu(square(newX.row(r))-square(p.X[l].row(r))) ;
        
            if( runif(1)[0] < exp( logLikRatio ) ){
                p.Xaccept[l](r) +=1;
            
                p.theta[l]  = newTheta;
                p.X[l]      = newX;
                p.logLik[l] = newLogLik;
            }
        }
    }
}





