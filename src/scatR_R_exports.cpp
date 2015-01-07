#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#include "scatR_cov.h"
#include "scatR_util.h"

RcppExport SEXP read_allele_file(SEXP Rfile) {

    std::string file_name = Rcpp::as<std::string>(Rfile);

    std::ifstream filestream(file_name.c_str(), std::ios_base::in | std::ios_base::binary);

    if (filestream.is_open()) {
        std::vector<char> data;
        try {
            boost::iostreams::filtering_istream f;
            f.push( boost::iostreams::gzip_decompressor() );
            f.push( filestream );

            boost::iostreams::copy(f, boost::iostreams::back_inserter(data));
        } catch (const boost::iostreams::gzip_error& e) {
            filestream.seekg(std::ios_base::beg);
            boost::iostreams::copy(filestream, boost::iostreams::back_inserter(data));
        }
        arma::vec x( (double *) &data[0], data.size() / sizeof(double), false);

        return(Rcpp::wrap(x));

    } else {
        Rcpp::Rcout << "Unable to open allele file (" << file_name <<")" << std::endl;
        return(R_NilValue);
    }
}

RcppExport SEXP R_cov_matern( SEXP rsigma2, SEXP rphi, SEXP rnu, SEXP rnugget, SEXP rdist, SEXP rdistmat, SEXP ruselog)
{
    double nugget = Rcpp::as<double>(rnugget);
    double sigma2 = Rcpp::as<double>(rsigma2);
    double nu     = Rcpp::as<double>(rnu);
    double phi    = Rcpp::as<double>(rphi);
    bool distmat  = Rcpp::as<bool>(rdistmat);
    bool uselog   = Rcpp::as<bool>(ruselog);

    SEXP res;
    if (distmat) {
        arma::mat dist = Rcpp::as<arma::mat>(rdist);

        res = Rcpp::wrap(cov_matern(sigma2, nu, phi, nugget, dist, uselog));
    } else {
        //NumericVector tDist(rdist);
        arma::vec dist = Rcpp::as<arma::vec>(rdist);

        res = Rcpp::wrap(cov_matern_vec(sigma2, nu, phi, nugget, dist, uselog));
    }

    return( res );
}

RcppExport SEXP R_great_circle_dist(SEXP x, SEXP y)
{
    return Rcpp::wrap( distance_mat(Rcpp::as<arma::colvec>(x), Rcpp::as<arma::colvec>(y)) );
}

RcppExport SEXP R_great_circle_dist_point(SEXP px, SEXP py, SEXP x, SEXP y)
{
    Rcpp::NumericVector xc(x), yc(y);
    Rcpp::NumericVector dist(xc.size());

    for(int i=0; i < xc.size(); ++i) {
        dist[i] = great_circle_dist(Rcpp::as<double>(px), Rcpp::as<double>(py), xc[i], yc[i]);
    }

    return(dist);
}

RcppExport SEXP R_fast_aggregate(SEXP R_rast, SEXP R_rstep, SEXP R_cstep, SEXP R_nrow, SEXP R_ncol, SEXP R_nlayer)
{
    int rstep  = Rcpp::as<int>(R_rstep),
        cstep  = Rcpp::as<int>(R_cstep),
        ncol   = Rcpp::as<int>(R_ncol),
        nrow   = Rcpp::as<int>(R_nrow),
        nlayer = Rcpp::as<int>(R_nlayer);

    int outrow = nrow/rstep,
        outcol = ncol/cstep;

    Rcpp::NumericMatrix rast(R_rast),
                        out( outrow*outcol, nlayer);


    for(int l=0; l < nlayer; l++) {

        for(int r=0; r<outrow; r++) {
            for(int c=0; c<outcol; c++) {

                double sum = 0.0;
                int n = 0;
                for(int rs=r*rstep; rs<(r+1)*rstep; rs++) {
                    for(int cs=c*cstep; cs<(c+1)*cstep; cs++) {
                        double val = rast(rs+cs*nrow,l);
                        if ( !R_IsNA(val) ) {
                            sum += val;
                            n++;
                        }
                    }
                }
                out(r+c*outrow, l) = sum / n;
            }
        }
    }

    return Rcpp::wrap( out );
}