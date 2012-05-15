#include <RcppArmadillo.h>
#include <boost/lexical_cast.hpp>

#include "Rscat_structs.h"
#include "Rscat_init.h"
#include "Rscat_mcmc.h"
#include "Rscat_util.h"
#include "Rscat_locate.h"

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


        //p.alfileStreams.resize(p.nLoci);
        //for(int l=0; l<p.locate_indivs.size(); l++) {
        //
        //    std::string file = opt.TMPDIR + "/CV_Ind" 
        //                     + boost::lexical_cast<std::string>(p.locate_indivs[l]) + "_"
        //                     + boost::lexical_cast<std::string>(p.chain_num) + ".gz";
        //    
        //    p.cvfileStreams.push_back( gzip_stream(file) );
        //    
        //    for(int j=0; j<p.nAlleles[l]; j++) {
        //        std::string al_file = opt.TMPDIR + "/Al" 
        //                            + boost::lexical_cast<std::string>(l+1) + "-"
        //                            + boost::lexical_cast<std::string>(j+1) + "_"
        //                            + boost::lexical_cast<std::string>(p.chain_num) + ".gz";
        //
        //        p.alfileStreams[l].push_back( gzip_stream(al_file) );
        //    }
        //}
        
        open_cvfiles(p,opt);
        open_allelefiles(p,opt);
        init_locate(p, opt);
    }
    
    Rcpp::List res = MCMCLoop(p, opt, Niter, Nthin, false, false);
    
    if (opt.VERBOSE) {
        std::cout << "Chain " << p.chain_num << ":" << std::endl;
        outputAccepts(p, opt);
        std::cout << std::endl << std::endl;
    }

    if (opt.LOCATE) {
          close_cvfiles(p, opt);
          close_allelefiles(p, opt);
    }
    
    return(res);
}