#include <RcppArmadillo.h>
#include <boost/lexical_cast.hpp>

#include "scatR_structs.h"
#include "scatR_init.h"
#include "scatR_mcmc.h"
#include "scatR_util.h"
#include "scatR_locate.h"

#ifdef USEMAGMA
#include "scatR_gpu.h"
#endif

RcppExport
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
    
    Rcpp::RNGScope scope;
    
    GlobalParams p;
    GlobalOptions opt;
    
    parseOptions(rOpt,opt);
    
    p.chain_num = Rcpp::as<int>(rChain);
    
    int Niter = Rcpp::as<int>(rNiter);
    int Nthin = Rcpp::as<int>(rNthin);
    int Nburn = Rcpp::as<int>(rNburn);
    
    p.indRegion = Rcpp::IntegerVector(rIndRegion);
    
    p.genotypes = Rcpp::as<arma::imat>(rGenotypes);
    p.boundary  = Rcpp::as<arma::mat>(rBoundary);
    p.locs      = Rcpp::as<arma::mat>(rLocations);
    
    
    p.nAlleles = Rcpp::IntegerVector(rNalleles);
    p.nRegions = p.locs.n_rows;
    p.nLoci = p.genotypes.n_cols;
    p.nInd = p.genotypes.n_rows/2;
    
    init_params(p,opt);
    
    std::cout << "Chain " << p.chain_num << ":\n";
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
        p.locate_indivs = Rcpp::IntegerVector(rLocIndivs);
        p.locate_genotypes = Rcpp::as<arma::imat>(rLocGenotypes);

        
        open_cvfiles(p,opt);
        if (opt.OUTPUTALFREQ)
            open_allelefiles(p,opt);
        init_locate(p, opt);
        
        #ifdef USEMAGMA
        init_GPU_data(p);
        #endif
    }
    
    Rcpp::List res = MCMCLoop(p, opt, Niter, Nthin, false, false);
    
    if (opt.VERBOSE) {
        std::cout << "Chain " << p.chain_num << ":" << std::endl;
        outputAccepts(p, opt);
        
        outputPerformance(p, opt);
    }

    if (opt.LOCATE) {
        close_cvfiles(p, opt);
        if (opt.OUTPUTALFREQ)
            close_allelefiles(p, opt);
        
        #ifdef USEMAGMA
        cleanup_GPU_data(p);
        #endif
    }
    
    return(res);
}