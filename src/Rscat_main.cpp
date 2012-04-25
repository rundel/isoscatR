#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

SEXP mcmc_main(SEXP rChain,
               SEXP rBoundary,   // px2 matrix
               SEXP rLocations,  // Rx2 matrix
               SEXP rGenotypes, // 2Ix(L+1) matrix
               SEXP rIndRegion,
               SEXP rNalleles,
               SEXP rNiter,
               SEXP rNthin,
               SEXP rNburn,
               SEXP rLocIndivs,
               SEXP rLocGenotypes,
               SEXP rOpt ) {
    
    RNGScope scope;
    
    GlobalParams p;
    GlobalOptions opt;
    
    parseOptions(rOpt,opt);
    
    p.chain_num = as<int>(rChain);
    
    int Niter = as<int>(rNiter);
    int Nthin = as<int>(rNthin);
    int Nburn = as<int>(rNburn);
    
    p.indRegion = IntegerVector(rIndRegion);
    
    
    IntegerMatrix tG(rGenotypes);  
    NumericMatrix tB(rBoundary);
    NumericMatrix tL(rLocations);
    p.genotypes = imat(tG.begin(),tG.nrow(),tG.ncol(),false);
    p.boundary = mat(tB.begin(),tB.nrow(),tB.ncol(),false);
    p.locs = mat(tL.begin(),tL.nrow(),tL.ncol(),false);
    
    
    p.nAlleles = IntegerVector(rNalleles);
    p.nRegions = p.locs.n_rows;
    p.nLoci = p.genotypes.n_cols;
    p.nInd = p.genotypes.n_rows/2;
    
    init_params(p,opt);
    
    
    cout << "Chain " << p.chain_num << ":\n";
    MCMCLoop(p, opt, Nburn, Nthin, true, true);
    if (opt.VERBOSE) {
        outputTuning(p, opt);
    }
    
    MCMCLoop(p, opt, Nburn, Nthin, true, false);    
    if (opt.VERBOSE) {
        outputAccepts(p, opt);
    }
    
    
    init_attempts(p);
    
	if (opt.LOCATE) {
        p.locate_indivs = IntegerVector(rLocIndivs);
        
        IntegerMatrix tmp(rLocGenotypes);  
        p.locate_genotypes = imat(tmp.begin(),tmp.nrow(),tmp.ncol(),false);   

        open_cvfiles(p, opt);
        open_allelefiles(p, opt);
        init_locate(p, opt);
    }
    
    List res = MCMCLoop(p, opt, Niter, Nthin, false, false);
    
    if (opt.VERBOSE) {
        cout << "Chain " << p.chain_num << ":" << endl;
        outputAccepts(p, opt);
        cout << endl << endl;
    }
    
    if (opt.LOCATE) {
        close_cvfiles(p, opt);
        close_allelefiles(p, opt);
    }
    
    return(res);
}