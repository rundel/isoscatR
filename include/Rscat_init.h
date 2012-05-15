#ifndef _RSCAT_INIT_H
#define _RSCAT_INIT_H

#include <RcppArmadillo.h>
#include "Rscat_structs.h"

void parseArgs(GlobalParams &p, SEXP rBoundary, SEXP rLocations, SEXP rRegionNames,
                                SEXP rGenotypes, SEXP rIndivID, SEXP rNalleles);
void parseOptions(SEXP sexpRopt, GlobalOptions &opt);

void init_params(GlobalParams &p, GlobalOptions &opt);
void init_proposal_sd(GlobalParams &p, GlobalOptions &opt);
void init_attempts(GlobalParams &p);
void calc_counts(GlobalParams &p, GlobalOptions &opt);
void calc_params(GlobalParams &p, GlobalOptions &opt);

#endif