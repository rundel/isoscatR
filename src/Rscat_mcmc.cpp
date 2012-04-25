#include <vector>
#include <boost/lexical_cast.hpp> 
#include "Rscat.h"
#include "Rscat_progressbar.h"

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

    mat samp_alpha, samp_mu, samp_xi, samp_beta, samp_aniso;
    vector<mat> samp_eta(p.nLoci);
    vector< vector<mat> > samp_X(p.nLoci);
    mat theta_fit;

    string prefix = tune ? "Tuning" : (burnin ? "Burnin" : "Sampling");

    progress_display progress_bar(Niter, cout, prefix);

    if(!burnin) {
        
        samp_alpha = mat(Niter,p.alpha.size());
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
		samp_aniso = mat(Niter,2);

        theta_fit = mat(Niter,1 + p.nRegions);
    }

    for (int i=0; i< Niter; i++){

        for (int j=0; j < Nthin; j++){
            MCMCStep(p, opt, burnin);
        }
        
        if(burnin && tune && i*Nthin % opt.TUNEINTERVAL == 0) 
            MCMCTune(p,opt);
        
        if (!burnin) {
            
            if (opt.LOCATE)
                update_Location(p, opt);
            
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
                for(int a=0; a<p.nAlleles[l]; a++) {
                    samp_X[l][a].row(i) = trans(p.X[l].col(a));
                }
                
                //deviance.row(i) += -2*trans(p.logLik[l]);
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
        
        ++progress_bar;
    }
    
    List res;
    
    if (!burnin) {
        
        vector<string> names_alpha, names_aniso, names_xibeta;
        names_alpha.push_back("alpha[0] - sigma2");
        names_alpha.push_back("alpha[1] - phi");
        names_alpha.push_back("alpha[2] - nu");
        names_alpha.push_back("alpha[3] - tau");
        
        names_aniso.push_back("anisotropy angle");
        names_aniso.push_back("anisotropy ratio");
        
        for (int l=0; l < p.nLoci; l++)
            names_xibeta.push_back("|xi|*beta [" + boost::lexical_cast<string>(l) + "]");
        
        
        res["alpha"] = List::create(Named("names") = names_alpha,
                                    Named("values") = samp_alpha);
                                    
        res["xibeta"] = List::create(Named("names") = names_xibeta,
                                     Named("values") = abs(samp_xi) % samp_beta);
        
        //res["deviance"] = List::create(Named("names") = "deviance",
        //                               Named("values") = sum(deviance,1));
        
        if (!opt.FIXRATIO && !opt.FIXANGLE) {
            res["aniso"] = List::create(Named("names") = names_aniso,
                                        Named("values") = samp_aniso);
        }
                                 
        if (opt.RETURNFIT) {
            
            vector<string> names_theta;
            names_theta.push_back("theta rss - total");
            
            for(int r=0; r<p.nRegions; ++r)
                names_theta.push_back("theta rss - Reg " + boost::lexical_cast<string>(r) ); 
            
            
            res["theta"] = List::create(Named("names") = names_theta,
                                        Named("values") = theta_fit);
        }
        
        // Calculating DIC based on fomula 6.12 on pg 183 of Bayesian Data Analysis
        //
        //vector<double> mean_alpha(p.alpha.size());
        //for(int i=0; i<p.alpha.size(); i++) {
        //    mean_alpha[i] = mean(samp_alpha.col(i));
        //}
        //
        //mat mean_S = calc_Sigma(mean_alpha, p.dist, opt.USEMATERN);
        //mat mean_L = calc_L(mean_S);
        //
        //colvec D_theta_hat = zeros<colvec>(p.nRegions);
        //for (int l=0; l < p.nLoci; l++) {
        //    double mean_xi = mean(samp_xi.col(l));
        //    
        //    rowvec mean_eta = mean(samp_eta[l],0);
        //    mat mean_X(p.nRegions,p.nAlleles[l]);
        //    
        //    for(int a=0; a<p.nAlleles[l]; a++) {
        //        mean_X.col(a) = trans(mean(samp_X[l][a],0));
        //    }        
        //    mat mean_theta = ones<colvec>(p.nRegions) * (mean_xi * mean_eta) + mean_L * mean_X;
        //
        //    D_theta_hat += -2*calc_multinomial_loglik(mean_theta, p.count[l], p.sumCount[l]);
        //}
        //
        //
        //colvec D_hat = trans(mean(deviance,0));
        //colvec pV = trans(var(deviance,0,0)/2);
        //colvec pD = D_hat - D_theta_hat;
        //
        //cout << "D_hat by region: ";
        //cout << trans(D_hat);
        //cout << "D_theta_hat by region: ";
        //cout << trans(D_theta_hat);
        
        //if (opt.VERBOSE) {
        //    cout << "Deviance Results:" << endl;
        //    cout << "=============================================" << endl;
        //
        //    cout << "D hat: " << setprecision(6) << setw(8) << accu(D_hat) << endl;
        //    cout << "pD: " << setprecision(6) << setw(8) << accu(pD) << endl;
        //    cout << "pV: " << setprecision(6) << setw(8) << accu(pV) << endl;
        //
        //    cout << "DIC total (pD): " << setprecision(6) << setw(8) << accu(D_hat+pD) << endl;
        //    cout << "DIC by region:";
        //    for(int r=0; r<p.nRegions; r++)
        //        cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pD(r)+0.5);
        //    cout << endl;
        //
        //    cout << "DIC total (pV): " << setprecision(6) << setw(8) << accu(D_hat+pV) << endl;
        //    cout << "DIC by region:";
        //    for(int r=0; r<p.nRegions; r++)
        //        cout << " " << setprecision(3) << setw(5) << floor(D_hat(r)+pV(r)+0.5);
        //    cout << endl << endl;
        //}
        
    }
    return(res);
}


void MCMCStep( GlobalParams &p, GlobalOptions &opt, bool burnin) {
    
    update_Beta(p, opt);    //cout << "beta" << endl;
    update_Eta(p,opt);      //cout << "eta" << endl;
    update_Xi(p,opt);       //cout << "xi" << endl;
    update_Alpha(p, opt);   //cout << "alpha" << endl;
    update_X(p, opt);       //cout << "X" << endl;

    //update_anisotropy(p,opt);
}

/*
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
            newTheta[l] = calc_theta(p.eta[l], p.xi[l], newL, p.X[l]);
            
            newLogLik[l] = calc_multinomial_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += accu(newLogLik[l] - p.logLik[l]);
        }
        
        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
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
        double newAngle = fmod(p.anisoAngle + rnorm(1,0,p.angle_sd)[0],  2*math::pi());

        mat newAngleMat = calc_rotation_mat (newAngle);
        mat newLocs = p.locs * newAngleMat * p.anisoRatioMat;
        mat newDist = calc_distance_mat(newLocs);
        mat newL    = calc_L(p.alpha, newDist, opt.USEMATERN);

        vector<mat> newTheta(p.nLoci);
        vector<colvec> newLogLik(p.nLoci);

        double logLikRatio = 0;
        for(int l=0; l<p.nLoci; l++) {
            newTheta[l] = calc_theta(p.eta[l], p.xi[l], newL, p.X[l]);

            newLogLik[l] = calc_multinomial_loglik(newTheta[l], p.count[l], p.sumCount[l]);

            logLikRatio += accu(newLogLik[l] - p.logLik[l]);
        }

        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
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


void update_Alpha(GlobalParams &p, GlobalOptions &opt) {   // 4              
    

    for(int i = 0; i < 4; i++){
        if (opt.FIXALPHA[i]) continue;
        
        p.alphaAttempt[i] +=1;
        
        vector<double> newAlpha(p.alpha);
        mat newL, newSigma;
        do {
            double h = rnorm(1,0,p.alpha_sd(i))[0];
            
            if (i == 1)
                newAlpha[i] = p.alpha[i] * exp( h );
            else
                newAlpha[i] = p.alpha[i] + h;
            
            
            if (newAlpha[i] < opt.ALPHAMIN[i] || newAlpha[i] > opt.ALPHAMAX[i])
                continue;
            
            newSigma = calc_Sigma(newAlpha, p.dist, opt.USEMATERN);
            newL = calc_L(newSigma);
            
        } while (newL.is_empty());
        
        mat newLInv = inv( trimatl(newL) );
        mat newSigmaInv = newLInv.t() * newLInv;
        double newSigmaDet = det(newSigma);
        
        // FIXME - transition probabilities are not symmetric 
        
        vector<double> newLogLik_theta(p.nLoci);
        
        double logLikRatio = 0;
        
        
        if (i == 0) {
            logLikRatio += -log(newAlpha[i]/p.alpha[i])-1/newAlpha[i]+1/p.alpha[i];
        }
        if (i == 3) {
            logLikRatio += -log(newAlpha[i]/p.alpha[i])-0.2*(1/newAlpha[i]-1/p.alpha[i]);
        }
        
		//cout << setprecision(5) << setw(7);
        //cout << "\nalpha[" << i << "]: " << p.alpha[i] << " -> " << newAlpha[i] << "\n";
        //cout << "sigmadet: " << newSigmaDet << " " << log(newSigmaDet) << "\n";
        
        for(int l=0; l<p.nLoci; l++) {
            
            mat mean = ones<colvec>(p.X[l].n_rows) * (p.xi[l] * p.eta[l]);
			newLogLik_theta[l] = calc_multivar_normal_loglik(p.theta[l], mean, newSigmaInv, newSigmaDet);            

			logLikRatio += newLogLik_theta[l]-p.logLik_theta[l];

            //cout << "theta: " << newLogLik_theta << " " << p.logLik_theta[l] << " ("
            //                  << calc_multivar_normal_loglik(p.theta[l], mean, p.Sinv, p.Sdet) << ") = "
			//				  << newLogLik_theta-p.logLik_theta[l] << "\n";
            
           
            
        }
        
        //cout << "llr: " << logLikRatio << " " << exp(logLikRatio) << "\n\n";

        
        if( runif(1)[0] < exp(logLikRatio) ) { //accept move 
            p.alphaAccept[i] +=1;
            
            p.alpha[i]     = newAlpha[i];
            p.L            = newL;
            p.Sinv         = newSigmaInv;
            p.Sdet         = newSigmaDet;
            //p.theta        = newTheta;
            p.logLik_theta = newLogLik_theta;
            //p.logLik_f     = newLogLik_f;
        }
        
    }
}

void update_X(GlobalParams &p, GlobalOptions &opt) {
    
    // Based on Eqn A.2 from BDA
    
    for(int l=0; l<p.nLoci; l++) {
        
        mat mean = ones<colvec>(p.nRegions) * (p.xi[l] * p.eta[l]);
        mat norm = randn<mat>(p.nRegions,p.nAlleles[l]);
        
        for(int r=0; r<p.nRegions; r++) {
            
            mat newTheta = p.theta[l];
            
            newTheta.row(r) = p.theta[l].row(r) + p.x_sd[l](r)*norm.row(r);
            
            p.Xattempt[l](r) += 1;
            
            double newLogLik_theta = calc_multivar_normal_loglik(newTheta, mean, p.Sinv, p.Sdet);
            double newLogLik_f     = calc_multinomial_loglik(newTheta, p.count[l], p.sumCount[l]);
            
            double logLikRatio = newLogLik_theta-p.logLik_theta[l]+newLogLik_f-p.logLik_f[l];
            
            if(0 && r==0 && l==0) {
                cout << "f:     " << newLogLik_f << " " << p.logLik_f[l] << " ("
                                  << calc_multinomial_loglik(p.theta[l], p.count[l], p.sumCount[l]) << ") = "
                                  << newLogLik_f-p.logLik_f[l] << "\n";            
                cout << "theta: " << newLogLik_theta << " " << p.logLik_theta[l] << " ("
                                  << calc_multivar_normal_loglik(p.theta[l], mean, p.Sinv, p.Sdet) << ") = "
                                  << newLogLik_theta-p.logLik_theta[l] << "\n"; 
                cout << "llr: " << logLikRatio << " (" << exp(logLikRatio) << ")\n\n";
            }
            if( runif(1)[0] < exp( logLikRatio ) ){
                p.Xaccept[l](r) +=1;
            
                p.theta[l]        = newTheta;
                //p.X[l]            = newX;
                p.logLik_theta[l] = newLogLik_theta;
                p.logLik_f[l]     = newLogLik_f;
            }
            
        }
        
        
        //mat newX = p.X[l] + p.x_sd[l](r) * randn<mat>(p.X[l].n_rows,p.X[l].n_cols);
        
        

    
        
        /*for(int r=0; r<p.nRegions; r++) {
            p.Xattempt[l](r) += 1;
        
            mat newX = p.X[l];
            newX.row(r) += p.x_sd[l](r) * randn<rowvec>(p.nAlleles[l]);
            
            mat newTheta = p.L * newX + mean;
            
            double newLogLik_theta = calc_multivar_normal_loglik(newTheta, mean, p.Sinv, p.Sdet);
            double newLogLik_f     = calc_multinomial_loglik(newTheta, p.count[l], p.sumCount[l]);
            
            double logLikRatio = newLogLik_theta-p.logLik_theta[l]+newLogLik_f-p.logLik_f[l]
                                 - 0.5 * accu( newX.row(r)%newX.row(r)-p.X[l].row(r)%p.X[l].row(r) );
            if(r==0 && l==0) {
            cout << "f:     " << newLogLik_f << " " << p.logLik_f[l] << " ("
                              << calc_multinomial_loglik(p.theta[l], p.count[l], p.sumCount[l]) << ") = "
                              << newLogLik_f-p.logLik_f[l] << "\n";            
            cout << "theta: " << newLogLik_theta << " " << p.logLik_theta[l] << " ("
                              << calc_multivar_normal_loglik(p.theta[l], mean, p.Sinv, p.Sdet) << ") = "
                              << newLogLik_theta-p.logLik_theta[l] << "\n"; 
			cout << "X:     " << - 0.5 * accu( newX.row(r)%newX.row(r)-p.X[l].row(r)%p.X[l].row(r) ) << "\n";
            cout << "llr: " << logLikRatio << " (" << exp(logLikRatio) << ")\n\n";
            }
            if( runif(1)[0] < exp( logLikRatio ) ){
                p.Xaccept[l](r) +=1;
                
                p.theta[l]        = newTheta;
                p.X[l]            = newX;
                p.logLik_theta[l] = newLogLik_theta;
                p.logLik_f[l]     = newLogLik_f;
            }
        }*/
    }
}


void update_Xi(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXXI)
        return;
    
    for(int l=0; l < p.nLoci; l++) {
        
        p.xiAttempt[l] +=1;
        
        double newXi = p.xi[l] + rnorm(1,0,p.xi_sd(l))[0];
        
        //mat old_mean = ones<colvec>(p.X[l].n_rows) * (p.xi[l] * p.eta[l]);
        mat new_mean = ones<colvec>(p.X[l].n_rows) * (newXi   * p.eta[l]);
        
		double newLogLik_theta = calc_multivar_normal_loglik(p.theta[l], new_mean, p.Sinv, p.Sdet);
		
        double logLikRatio = newLogLik_theta-p.logLik_theta[l];
        
        if( runif(1)[0] < exp(logLikRatio) ){ //accept move 
            p.xiAccept[l] += 1;

            p.xi[l]           = newXi;                
            p.logLik_theta[l] = newLogLik_theta;
        }
    }
}



void update_Eta(GlobalParams &p, GlobalOptions &opt) {    // NA  1xA[l]
    
    if (opt.FIXETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        p.etaAttempt[l] += 1;
        
        rowvec newEta = p.eta[l] + p.eta_sd(l) * randn<rowvec>(p.nAlleles[l]);
        
        //mat old_mean = ones<colvec>(p.X[l].n_rows) * (p.xi[l] * p.eta[l]);
        mat new_mean = ones<colvec>(p.X[l].n_rows) * (p.xi[l] * newEta  );
        
        double newLogLik_theta = calc_multivar_normal_loglik(p.theta[l], new_mean, p.Sinv, p.Sdet);
		
        //prior on eta is N(0,beta) 
        double logLikRatio = newLogLik_theta - p.logLik_theta[l] - 0.5*accu(square(newEta) - square(p.eta[l])) / (p.beta[l]*p.beta[l]);
        
        if( runif(1)[0] < exp(logLikRatio) ){ //accept move 
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

void update_Beta(GlobalParams &p, GlobalOptions &opt) {
    
    if (opt.FIXBETA) return;
    
    for(int l=0; l<p.nLoci; l++) {
        
        p.betaAttempt[l] += 1;
        
        double newBeta;
        do {
            newBeta = p.beta[l] + rnorm(1,0,p.beta_sd(l))[0];
        } while(newBeta < opt.BETARANGE[0] || newBeta > opt.BETARANGE[1]);
        
        double etasq_sum = accu(p.eta[l] % p.eta[l]);
        
        double logLikRatio = -p.nAlleles[l]*(log(newBeta/p.beta[l])) - 0.5 * etasq_sum * (1/(newBeta*newBeta) - 1/(p.beta[l]*p.beta[l]));
        
        if( runif(1)[0] < exp( logLikRatio ) ){ //accept move 
            p.betaAccept[l] += 1;
            
            p.beta[l] = newBeta;
        }
    }
}
