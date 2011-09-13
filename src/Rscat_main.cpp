#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;




SEXP mcmc_main(SEXP rChain,
              SEXP rBoundary,   // px2 matrix
              SEXP rLocations,  // Rx2 matrix
              SEXP rGeneotypes, // 2Ix(L+1) matrix
              SEXP rIndRegion,
              SEXP rNalleles,
              SEXP rNiter,
              SEXP rNthin,
              SEXP rNburn,
              SEXP rCVIndivs,
              SEXP rCVGeneotypes,
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
    
    
    IntegerMatrix tG(rGeneotypes);  
    NumericMatrix tB(rBoundary);
    NumericMatrix tL(rLocations);
    p.geneotypes = imat(tG.begin(),tG.nrow(),tG.ncol(),false);
    p.boundary = mat(tB.begin(),tB.nrow(),tB.ncol(),false);
    p.locs = mat(tL.begin(),tL.nrow(),tL.ncol(),false);
    
    
    p.nAlleles = IntegerVector(rNalleles);
    p.nRegions = p.locs.n_rows;
    p.nLoci = p.geneotypes.n_cols;
    p.nInd = p.geneotypes.n_rows/2;
    
    init_params(p,opt);
    
    MCMCLoop(p, opt, Nburn, Nthin, true, true);
    if (opt.VERBOSE)
        outputTuning(p, opt);
    
    MCMCLoop(p, opt, Nburn, Nthin, true, false);    
    
    init_attempts(p);
    if (opt.LOCATE) {
        init_locate(p, opt);
        open_allelefiles(p, opt);
    }
    if (opt.CROSSVALIDATE) {
        p.cvIndivs = IntegerVector(rCVIndivs);
        
        IntegerMatrix tCVG(rCVGeneotypes);  
        p.cvGeneotypes = imat(tCVG.begin(),tCVG.nrow(),tCVG.ncol(),false);
        
        
        init_locate(p, opt);
        open_cvfiles(p, opt);
    }
    
    List res = MCMCLoop(p, opt, Niter, Nthin, false, false);
    
    if (opt.VERBOSE) {
        cout << "Chain " << p.chain_num << ":" << endl;
        outputTuning(p, opt);
        outputAccepts(p, opt);
        cout << endl << endl;
    }
    
    if (opt.LOCATE) close_allelefiles(p, opt);
    if (opt.CROSSVALIDATE) close_cvfiles(p, opt);
    
    return(res);
}