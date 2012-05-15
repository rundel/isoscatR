#ifndef _RSCAT_MCMC_H
#define _RSCAT_MCMC_H

#include "Rscat_structs.h"

double tuneScale(double acc_rate);
void MCMCTune(GlobalParams &p, GlobalOptions &opt);

Rcpp::List MCMCLoop(GlobalParams &p, GlobalOptions &opt, int Niter, int Nthin, bool burnin, bool thin);
void MCMCStep( GlobalParams &p, GlobalOptions &opt, bool burnin);    
void update_beta(GlobalParams &p, GlobalOptions &opt);
void update_theta(GlobalParams &p, GlobalOptions &opt);
void update_alpha(GlobalParams &p, GlobalOptions &opt);
void update_eta(GlobalParams &p, GlobalOptions &opt);
void update_xi(GlobalParams &p, GlobalOptions &opt);
void update_anisotropy(GlobalParams &p, GlobalOptions &opt);

#endif