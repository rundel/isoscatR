#include <iostream>
#include <RcppArmadillo.h>

#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>

#include "scatR_structs.h"
#include "scatR_util.h"
#include "scatR_cov.h"

#ifdef USEMAGMA
#include "scatR_magma_wrap.h"
#include "scatR_gpu.h"
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

    //if (opt.VERBOSE) {
        Rcpp::Rcout << "Estimated # grid cells: " << nx_step*ny_step << std::endl;
        Rcpp::Rcout << "Actual # grid cells: " << i-1 << std::endl;
        Rcpp::Rcout << "Using step=" << step << std::endl;
        Rcpp::Rcout << std::endl;
    //}

    arma::mat rand_pred_locs = arma::shuffle(pred_locs,0);
    arma::mat all_locs = arma::join_cols(p.locs, rand_pred_locs);

    p.pred_dist = distance_mat(all_locs);

    std::string file = opt.TMPDIR + "/pred_coords_"
                       + boost::lexical_cast<std::string>(p.chain_num) + ".mat";

    rand_pred_locs.save(file, arma::raw_ascii);

    arma::rowvec x = arma::trans(rand_pred_locs.col(0)),
                 y = arma::trans(rand_pred_locs.col(1));
}

void update_location(GlobalParams &p, GlobalOptions &opt) {

    boost::timer::cpu_timer timer;

    int n_known = p.nRegions;
    int n_total = p.pred_dist.n_rows;
    int n_pred = n_total-n_known;
    int total_alleles = Rcpp::sum(p.nAlleles);

#ifndef USEMAGMA

    timer.stop();
    timer.start();

    arma::mat dist12 = p.pred_dist(arma::span(0,n_known-1),        arma::span(n_known, n_total-1));
    arma::mat dist22 = p.pred_dist(arma::span(n_known, n_total-1), arma::span(n_known, n_total-1));

    std::vector<double> cur_alpha = p.alpha;
    cur_alpha[3] = 0;
    arma::mat cov12 = calc_Sigma(cur_alpha, dist12, opt.USEMATERN);
    arma::mat cov22 = calc_Sigma(p.alpha, dist22, opt.USEMATERN);

    p.step1(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    arma::mat tmp = cov12.t() * p.Sinv;

    p.step2(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    arma::mat predL = arma::chol(cov22 - tmp * cov12).t();

    p.step3(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    arma::mat pred = predL * arma::randn<arma::mat>(n_pred,total_alleles);

    p.step4(timer.elapsed().wall / 1000000000.0L);

#else

    double one = 1.0, negone = -1.0, zero = 0.0;


    timer.stop();
    timer.start();

    cov_powered_exponential_gpu(p.d_dist12, p.d_cov12, n_known, n_pred,  p.alpha[0], p.alpha[1], p.alpha[2], 0.0, 64);
    cov_powered_exponential_gpu(p.d_dist22, p.d_cov22, n_pred,  n_pred,  p.alpha[0], p.alpha[1], p.alpha[2], p.alpha[3], 64);

    checkCublasError(
        cublasSetMatrix(n_known, n_known, sizeof(double), p.Sinv.memptr(), n_known, p.d_invcov11, n_known),
        "Sinv (Set)"
    );

    p.step1(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    checkCublasError( // tmp = t(cov12) * cov11^-1
        cublasDgemm_v2( p.handle, CUBLAS_OP_T, CUBLAS_OP_N,
                        n_pred, n_known, n_known,
                       &one,
                        p.d_cov12, n_known,
                        p.d_invcov11, n_known,
                       &zero,
                        p.d_tmp, n_pred ),
        "tmp"
    );

    p.step2(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    arma::mat tmp(n_pred, n_known);
    checkCublasError( cublasGetMatrix(n_pred, n_known, sizeof(double), p.d_tmp, n_pred, tmp.memptr(), n_pred), "tmp (Get)" );

    checkCublasError( // cov22 = cov22 - tmp * cov12
        cublasDgemm_v2( p.handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n_pred, n_pred, n_known,
                       &negone,
                        p.d_tmp, n_pred,
                        p.d_cov12, n_known,
                       &one,
                        p.d_cov22, n_pred ),
        "predL"
    );

    int info;
    magma_cholesky('L', n_pred, p.d_cov22, n_pred, &info);

    if (info != 0) {

        if (info < 0) Rcpp::Rcout << "Error: illegal argument " << -1 * info << " in magChol\n";
        else          Rcpp::Rcout << "Error: leading minor of order " << info << " is not positive definite\n";

        Rcpp::Rcout << "alpha[0]: " << p.alpha[0] << "\n"
                  << "alpha[1]: " << p.alpha[1] << "\n"
                  << "alpha[2]: " << p.alpha[2] << "\n"
                  << "alpha[3]: " << p.alpha[3] << "\n";

        throw( std::exception() );
    }

    for(int col=1; col < n_pred; ++col) {
        double* colptr = p.d_cov22 + col*n_pred;
        cudaMemset(colptr, 0, col * sizeof(double));
    }

    p.step3(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    double *d_rand, *d_pred;
    cudaMalloc((void**) &d_rand, n_pred*total_alleles*sizeof(double));  checkCudaError("rand (Malloc)");
    cudaMalloc((void**) &d_pred, n_pred*total_alleles*sizeof(double));  checkCudaError("pred (Malloc)");

    arma::mat rand = arma::randn<arma::mat>(n_pred,total_alleles);
    checkCublasError(
        cublasSetMatrix(n_pred, total_alleles, sizeof(double), rand.memptr(), n_pred, d_rand, n_pred),
        "dist12 (Set)"
    );


    checkCublasError( // pred = predL * rand
        cublasDgemm_v2( p.handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n_pred, total_alleles, n_pred,
                       &one,
                        p.d_cov22, n_pred,
                        d_rand, n_pred,
                       &zero,
                        d_pred, n_pred ),
        "predL"
    );

    arma::mat pred(n_pred, total_alleles);
    checkCublasError(
        cublasGetMatrix(n_pred, total_alleles, sizeof(double), d_pred, n_pred, pred.memptr(), n_pred),
        "predL (Get)"
    );

    cudaFree(d_rand);  checkCudaError("rand (Free)");
    cudaFree(d_pred);  checkCudaError("pred (Free)");

    p.step4(timer.elapsed().wall / 1000000000.0L);

#endif

    timer.stop();
    timer.start();

    arma::mat ind_logprob = arma::zeros<arma::mat>(p.locate_indivs.size(), n_pred);
    for (int l=0, c=0; l<p.nLoci; c+=p.nAlleles[l], l++) {

        arma::mat mean_rows = p.xi[l] * p.eta[l];
        arma::mat mean      = (arma::ones<arma::colvec>(n_known) * mean_rows);
        arma::mat mean_pred = (arma::ones<arma::colvec>(n_pred)  * mean_rows) + tmp * (p.theta[l] - mean);

        //arma::mat rand = arma::randn<arma::mat>(n_pred,p.nAlleles[l]);
        //arma::mat res = mean_pred + predL * rand;
        arma::mat res = mean_pred + pred(arma::span::all, arma::span(c,c+p.nAlleles[l]-1));


        arma::mat alleleProb = (1-opt.DELTA) * calc_f( res ) + opt.DELTA/p.nAlleles[l];

        if (opt.OUTPUTALFREQ) {
            for(int j=0; j<p.nAlleles[l]; j++) {
                arma::rowvec curvec = trans(alleleProb.col(j));
                if (opt.GZIPOUTPUT) curvec.save( *(p.alfileGzStreams[l][j]), arma::raw_binary );
                else                curvec.save( *(p.alfileStreams[l][j]),   arma::raw_binary );
            }
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

    p.step5(timer.elapsed().wall / 1000000000.0L);


    timer.stop();
    timer.start();

    for(int i=0; i < p.locate_indivs.size(); i++) {
        arma::rowvec curvec = ind_logprob.row(i);
        if (opt.GZIPOUTPUT) curvec.save( *(p.cvfileGzStreams[i]), arma::raw_binary );
        else                curvec.save( *(p.cvfileStreams[i]),   arma::raw_binary );
    }

    p.step6(timer.elapsed().wall / 1000000000.0L);

}
