#include <iostream>
#include <RcppArmadillo.h>

#include <boost/lexical_cast.hpp> 

#include "Rscat_structs.h"
#include "Rscat_util.h"
#include "Rscat_cov.h"

#ifdef USEMAGMA
#include "Rscat_gpu.h"
#include <magma.h>
#endif


void init_locate(GlobalParams &p, GlobalOptions &opt) {

    double xmin = floor( min(p.boundary.col(0)) * 180.0/arma::math::pi());
    double ymin = floor( min(p.boundary.col(1)) * 180.0/arma::math::pi());
    double xmax =  ceil( max(p.boundary.col(0)) * 180.0/arma::math::pi());
    double ymax =  ceil( max(p.boundary.col(1)) * 180.0/arma::math::pi());

    int nx = (int) (xmax-xmin);
    int ny = (int) (ymax-ymin);
    
    double ratio = opt.MAXCELL / (nx*ny);
    double step;

    if (ratio > 1) {
        step = 1/floor(sqrt(ratio));
    } else {
        step = ceil(sqrt(1/ratio));
        
        nx += nx % (int) step;
        ny += ny % (int) step;
    }
    
    int nx_step = nx / step;
    int ny_step = ny / step;
    
    arma::mat pred_locs(nx_step*ny_step,2);
    
    int i = 0;
    for (int xi = 0; xi < nx_step; xi++) {
        for(int yi = 0; yi < ny_step; yi++) {
        
            double newX = (arma::math::pi()/180.0) * (xmin + (xi*step) + step/2.0);
            double newY = (arma::math::pi()/180.0) * (ymax - (yi*step) - step/2.0);
        
            if (isInsideBoundary(newX,newY, p.boundary) == 0)
                continue;
            
            pred_locs(i,0) = newX;
            pred_locs(i,1) = newY;
            
            i++;
        }
    }
    pred_locs = pred_locs(arma::span(0,i-1), arma::span::all); 
    
    if (opt.VERBOSE) {
        std::cout << "Estiarma::mated # grid cells: " << nx_step*ny_step << std::endl;
        std::cout << "Actual # grid cells: " << i-1 << std::endl;
        std::cout << "Using step=" << step << std::endl;
        std::cout << std::endl;
    }
    
    arma::mat rand_pred_locs = arma::shuffle(pred_locs,0);
    arma::mat all_locs = arma::join_cols(p.locs, rand_pred_locs);
    
    p.pred_dist = distance_mat(all_locs);
    
    std::string file = opt.TMPDIR + "/pred_coords_" 
                       + boost::lexical_cast<std::string>(p.chain_num) + ".mat";
    
    rand_pred_locs.save(file, arma::raw_ascii);
    
    arma::rowvec x = arma::trans(rand_pred_locs.col(0)),
                 y = arma::trans(rand_pred_locs.col(1));
    
    /*for (int l=0; l<p.nLoci; l++) {
        for(int j=0; j<p.nAlleles[l]; j++) {
            x.save( *(p.alfileGzStreams[l][j]), arma::raw_ascii );
            y.save( *(p.alfileGzStreams[l][j]), arma::raw_ascii );
        }
    }*/
    
#ifdef USEMAGMA
    p.gpu = GPU_data(p, opt);
#endif
}

void update_location(GlobalParams &p, GlobalOptions &opt) {
    
    arma::mat predL, tmp;
    
    int n_known = p.nRegions;
    int n_total = p.pred_dist.n_rows;
    int n_pred = n_total-n_known;
    
#ifndef USEMAGMA
    arma::mat dist12 = p.pred_dist(arma::span(0,n_known-1),        arma::span(n_known, n_total-1));
    arma::mat dist22 = p.pred_dist(arma::span(n_known, n_total-1), arma::span(n_known, n_total-1));
    
    arma::mat cov22 = calc_Sigma(p.alpha, dist22, opt.USEMATERN);
    
    std::vector<double> cur_alpha = p.alpha;
    cur_alpha[3] = 0;
    arma::mat cov12 = calc_Sigma(cur_alpha, dist12, opt.USEMATERN);
    
    tmp = cov12.t() * p.Sinv;
    predL = arma::chol(cov22 - tmp * cov12).t();
#else
    double one = 1.0, negone = -1.0, zero = 0.0;
    
    double *d_invcov11, *d_cov12, *d_cov22, *d_tmp;
    cudaMalloc((void**)&d_cov12, n_known*n_pred  * sizeof(double));
    cudaMalloc((void**)&d_cov22, n_pred *n_pred  * sizeof(double));
    cudaMalloc((void**)&d_invcov11, n_known*n_known * sizeof(double));
    cudaMalloc((void**)&d_tmp    , n_pred*n_known * sizeof(double));
    checkCudaError("Locate - Malloc");
    
    std::cout << "HERE1\n";
    //cov_powered_exponential_gpu(p.gpu.d_dist11, p.gpu.d_cov11, n_known, n_known, p.alpha[0], p.alpha[1], p.alpha[2], p.alpha[3], 64);
    cov_powered_exponential_gpu(p.gpu.d_dist12, d_cov12, n_known, n_pred,  p.alpha[0], p.alpha[1], p.alpha[2],        0.0, 64);
    cov_powered_exponential_gpu(p.gpu.d_dist22, d_cov22, n_pred,  n_pred,  p.alpha[0], p.alpha[1], p.alpha[2], p.alpha[3], 64);
    std::cout << "HERE2 " << n_known << " " <<  p.Sinv.n_rows << " " << p.Sinv.n_cols << "\n";
    checkCublasError( 
        cublasSetMatrix(n_known, n_known, sizeof(double), p.Sinv.memptr(), n_known, d_invcov11, n_known) 
        //cublasSetMatrix(n_known, n_known, sizeof(double), p.Sinv.memptr(), n_known, p.gpu.d_invcov11, n_known) 
    );
    std::cout << "HERE3\n";
    // tmp = t(cov12) * cov11^-1
    //checkCublasError( 
        cublasDgemm_v2( p.gpu.handle, CUBLAS_OP_T, CUBLAS_OP_N,
                        n_pred, n_known, n_known,
                       &one,
                        d_cov12, n_known,
                        d_invcov11, n_known, //p.gpu.d_invcov11, n_known,
                       &zero,
                        d_tmp, n_pred ) //p.gpu.d_tmp, n_pred )
    ;//);
    
    std::cout << "HERE4\n";
    
    // cov22 = cov22 - tmp * cov12
    //checkCublasError(
        cublasDgemm_v2( p.gpu.handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n_pred, n_known, n_pred,
                       &negone,
                        d_tmp, n_pred, //p.gpu.d_tmp, n_pred,
                        d_cov12, n_known,
                       &one,
                        d_cov22, n_pred )
    ;//);
    int info;
    magma_dpotrf_gpu('L', n_pred, d_cov22, n_pred, &info);
    
    if (info != 0) {
    
        if (info < 0)
            std::cout << "Error: illegal argument " << -1 * info << " in magChol\n"; 
        else
            std::cout << "Error: leading minor of order " << info << " is not positive definite\n";

        std::cout << "alpha[0]: " << p.alpha[0] << "\n"
                  << "alpha[1]: " << p.alpha[1] << "\n"
                  << "alpha[2]: " << p.alpha[2] << "\n"
                  << "alpha[3]: " << p.alpha[3] << "\n";

        throw( std::exception() );
    }
    
    std::cout << "HERE5\n";
    
    checkCublasError( cublasGetMatrix(n_pred, n_known, sizeof(double), d_tmp,   n_pred, tmp.memptr(),   n_pred) );
    checkCublasError( cublasGetMatrix(n_pred, n_pred,  sizeof(double), d_cov22, n_pred, predL.memptr(), n_pred) );
    
    for(int col=1; col < n_pred; ++col) {
        double* colptr = predL.colptr(col);
        memset(colptr, 0, col * sizeof(double));
    }
    
    cudaFree(d_cov12);
    cudaFree(d_cov22);
    cudaFree(d_invcov11);
    cudaFree(d_tmp);
    
    checkCudaError("Locate - Cleanup");
#endif
    
    arma::mat ind_logprob = arma::zeros<arma::mat>(p.locate_indivs.size(), n_pred);
    for (int l=0; l<p.nLoci; l++) {
        
        arma::mat mean_rows = p.xi[l] * p.eta[l];
        
        arma::mat mean      = (arma::ones<arma::colvec>(n_known) * mean_rows);
        arma::mat mean_pred = (arma::ones<arma::colvec>(n_pred)  * mean_rows) + tmp * (p.theta[l] - mean);
        
        arma::mat res = mean_pred + predL * arma::randn<arma::mat>(n_pred,p.nAlleles[l]);
        
        arma::mat alleleProb = (1-opt.DELTA) * calc_f( res ) + opt.DELTA/p.nAlleles[l];
        
        for(int j=0; j<p.nAlleles[l]; j++) {
            arma::rowvec curvec = trans(alleleProb.col(j));
            curvec.save( *(p.alfileStreams[l][j]), arma::raw_binary );
            //curvec.save( *(p.alfileStreams[l][j].stream()), arma::raw_binary );
        }
        
        for (int i=0; i < p.locate_indivs.size(); i++) {
            int al1 = p.locate_genotypes(2*i  , l);
            int al2 = p.locate_genotypes(2*i+1, l);

            arma::vec newP1, newP2;

            if (al1 < 0) newP1 = arma::ones<arma::vec>(n_pred);
            else         newP1 = alleleProb.col(al1);

            if (al2 < 0) newP2 = arma::ones<arma::vec>(n_pred);
            else         newP2 = alleleProb.col(al2);

            if (al1 == al2) ind_logprob.row(i) += arma::trans( log(opt.NULLPROB * newP1 + (1-opt.NULLPROB) * arma::square(newP1)) );
            else            ind_logprob.row(i) += arma::trans( log((1-opt.NULLPROB) * (newP1 % newP2)) );
        }
    }
    
    for(int i=0; i < p.locate_indivs.size(); i++) {
        arma::rowvec curvec = ind_logprob.row(i);
        curvec.save( *(p.cvfileStreams[i]), arma::raw_binary );
        //curvec.save( *(p.cvfileStreams[i].stream()), arma::raw_binary );
    }

}
