#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;


void outputAccepts(GlobalParams &p, GlobalOptions &opt) {
    
    cout << "Acceptance Rates:" << endl;
    cout << "=============================================" << endl;
    
    
    cout << "alpha  : ";
    for (int l=0; l<p.alpha.size(); l++) 
        cout << setprecision(3) << setw(5) << calc_accept_ratio(p.alphaAccept[l], p.alphaAttempt[l]) << " ";
    cout << endl;
    
    if (!opt.FIXANGLE) {
        cout << "angle  : ";
        cout << setprecision(3) << setw(5) << calc_accept_ratio(p.angleAccept, p.angleAttempt) << " ";
        cout << endl;
    }
    
    if (!opt.FIXRATIO) {
        cout << "ratio  : ";
        cout << setprecision(3) << setw(5) << calc_accept_ratio(p.ratioAccept, p.ratioAttempt) << " ";
        cout << endl;
    }
    
    if (!opt.FIXXI) {
        cout << "xi     : ";
        for (int l=0; l<p.nLoci; l++) 
            cout << setprecision(3) << setw(5) << calc_accept_ratio(p.xiAccept[l],p.xiAttempt[l]) << " ";
        cout << endl;
    }
    
    if (!opt.FIXBETA) {
        cout << "beta   : ";
        for (int l=0; l<p.nLoci; l++) 
            cout << setprecision(3) << setw(5) << calc_accept_ratio(p.betaAccept[l], p.betaAttempt[l]) << " ";
        cout << endl;
    }
    
    if (!opt.FIXETA) {
        cout << "eta    : ";
        for (int l=0; l<p.nLoci; l++) 
            cout << setprecision(3) << setw(5) << calc_accept_ratio(p.etaAccept[l], p.etaAttempt[l]) << " ";
        cout << endl;
    }
    
    for (int l=0; l<p.nLoci; l++) {
        cout << "X[" << l << "]   : ";
        for(int r=0; r<p.nRegions; r++) {
            cout << setprecision(3) << setw(5) << calc_accept_ratio(p.Xaccept[l](r), p.Xattempt[l](r)) << " ";
        }
        cout << endl;
    }
    cout << endl;
    
}


void outputTuning(GlobalParams &p, GlobalOptions &opt) {
    
    cout << endl;
    cout << "Tuning Results:" << endl;
    cout << "=============================================" << endl;

    
    for (int l=0; l<opt.ALPHASD.size(); l++) {
        cout << "alpha" << l << " : "; 
        cout << "[" << setprecision(3) << setw(5) << opt.ALPHASD(l) << "]";
        cout << " " << setprecision(3) << setw(5) << p.alpha_sd(l) << endl;
    }
    
    if (!opt.FIXANGLE) {
        cout << "angle  : ";
        cout << "[" << setprecision(3) << setw(5) << opt.ANGLESD << "]";
        cout << " " << setprecision(3) << setw(5) << p.angle_sd << endl;
    }
    
    if (!opt.FIXRATIO) {
        cout << "ratio  : ";
        cout << "[" << setprecision(3) << setw(5) << opt.RATIOSD << "]";
        cout << " " << setprecision(3) << setw(5) << p.ratio_sd << endl;
    }
    
    if (!opt.FIXXI) {
        cout << "xi     : [" << setprecision(3) << setw(5) << opt.XISD << "] ";
        cout << "(" << setprecision(3) << setw(5) << mean(p.xi_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            cout << " " << setprecision(3) << setw(5) << p.xi_sd(l);
        cout << endl;
    }
    
    if (!opt.FIXBETA) {
        cout << "beta   : [" << setprecision(3) << setw(5) << opt.BETASD << "] ";
        cout << "(" << setprecision(3) << setw(5) << mean(p.beta_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            cout << " " << setprecision(3) << setw(5) << p.beta_sd(l);
        cout << endl;
    }
    
    if (!opt.FIXETA) {
        cout << "eta    : [" << setprecision(3) << setw(5) << opt.ETASD << "] ";
        cout << "(" << setprecision(3) << setw(5) << mean(p.eta_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            cout << " " << setprecision(3) << setw(5) << p.eta_sd(l);
        cout << endl;
    }
    
    for (int l=0; l<p.nLoci; l++) {
        cout << "X[" << l << "]   : [" << setprecision(3) << setw(5) << opt.XSD << "] ";
        cout << "(" << setprecision(3) << setw(5) << mean(p.x_sd[l]) << ")";
        for(int r=0; r<p.nRegions; r++) {
            cout << " " << setprecision(3) << setw(5) << p.x_sd[l](r);
        }
        cout << endl;
    }
    cout << endl;
}






