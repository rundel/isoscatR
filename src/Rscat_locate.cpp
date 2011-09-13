#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

void init_locate(GlobalParams &p, GlobalOptions &opt) {

    double xmin = floor( min(p.boundary.col(0)) * 180.0/math::pi());
    double ymin = floor( min(p.boundary.col(1)) * 180.0/math::pi());
    double xmax =  ceil( max(p.boundary.col(0)) * 180.0/math::pi());
    double ymax =  ceil( max(p.boundary.col(1)) * 180.0/math::pi());

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
    
    mat pred_locs(nx_step*ny_step,2);
    
    //cout << endl;
    //cout << "nx=" << nx << endl;
    //cout << "ny=" << ny << endl;
    //cout << "nx_step=" << nx_step << endl;
    //cout << "ny_step=" << ny_step << endl;
    //cout << "xmn=" << xmin << ", xmx=" << xmax 
    //     << ", ymn=" << ymin << ", ymx=" << ymax << endl;
    
    int i = 0;
    for(int xi = 0; xi < nx_step; xi++) {
        for(int yi = 0; yi < ny_step; yi++) {
        
            double newX = (math::pi()/180.0) * (xmin + (xi*step) + step/2.0);
            double newY = (math::pi()/180.0) * (ymax - (yi*step) - step/2.0);
        
            if (isInsideBoundary(newX,newY, p.boundary) == 0)
                continue;
            
            pred_locs(i,0) = newX;
            pred_locs(i,1) = newY;
            
            i++;
        }
    }
    pred_locs = pred_locs(span(0,i-1), span::all); 
    
    cout << endl;
    cout << "Estimated # grid cells: " << nx_step*ny_step << endl;
    cout << "Actual # grid cells: " << i-1 << endl;
    cout << "Using step=" << step << endl;
    cout << endl;
    
    
    
    
    mat rand_pred_locs = shuffle(pred_locs,0);
    mat all_locs = join_cols(p.locs, rand_pred_locs);
    
    p.predDist = calc_distance_mat(all_locs);
    
    stringstream ss;
    ss << opt.TMPDIR << "/pred_coords";
    if (opt.CROSSVALIDATE)
        ss << p.cvIndivs[0] << opt.FILESUFFIX;
    ss << ".mat";

    rand_pred_locs.save(ss.str(),raw_ascii);
    
}

void update_Location(GlobalParams &p, GlobalOptions &opt) {
    
    
    int nPos = p.predDist.n_rows;
    int nknown = p.locs.n_rows;
    int npred = nPos - nknown;
    
    mat newL = calc_L(p.alpha,p.predDist,opt.USEMATERN);
 
    for(int l=0; l<p.nLoci; l++) {
        mat newX = join_cols(p.X[l],randn<mat>(npred,p.nAlleles[l]));
    
        mat mean = p.mu[l] * ones<mat>(npred,p.nAlleles[l]) + (ones<colvec>(npred) * (p.xi[l] * p.eta[l]));
        mat var = newL*newX;
        mat res = mean  + var(span(nknown,nPos-1), span::all);
        
        mat alleleProb = (1-opt.DELTA) * calc_f( res ) + opt.DELTA/p.nAlleles[l];
        
        for(int j=0; j<p.nAlleles[l]; j++) {
            rowvec curvec = trans(alleleProb.col(j));
            curvec.save( *(p.alfileGzStreams[l][j]), raw_binary );
        }
    }
}

void update_LocationCV(GlobalParams &p, GlobalOptions &opt) {
    
    
    int nPos = p.predDist.n_rows;
    int nknown = p.locs.n_rows;
    int npred = nPos - nknown;
    
    int nCVInd = p.cvIndivs.size();
    
    mat ind_logprob = zeros<mat>(nCVInd,npred);
    mat newL = calc_L(p.alpha,p.predDist,opt.USEMATERN);
 
    for(int l=0; l<p.nLoci; l++) {
        mat newX = join_cols(p.X[l],randn<mat>(npred,p.nAlleles[l]));
    
        mat mean = p.mu[l] * ones<mat>(npred,p.nAlleles[l]) + (ones<colvec>(npred) * (p.xi[l] * p.eta[l]));
        mat var = newL*newX;
        mat res = mean  + var(span(nknown,nPos-1), span::all);
        
        mat alleleProb = (1-opt.DELTA) * calc_f( res ) + opt.DELTA/p.nAlleles[l];
        
        for(int i=0; i < nCVInd; i++) {
            int al1 = p.cvGeneotypes(2*i  , l);
            int al2 = p.cvGeneotypes(2*i+1, l);

            vec newP1,newP2;

            if (al1 < 0) newP1 = ones<vec>(npred);
            else         newP1 = alleleProb.col(al1);

            if (al2 < 0) newP2 = ones<vec>(npred);
            else         newP2 = alleleProb.col(al2);

            if (al1 == al2) ind_logprob.row(i) += trans( log(opt.NULLPROB * newP1 + (1-opt.NULLPROB) * square(newP1)) );
            else            ind_logprob.row(i) += trans( log((1-opt.NULLPROB) * (newP1 % newP2)) );

        }
        
    }
    
    for(int i=0; i < nCVInd; i++) {
        rowvec curvec = ind_logprob.row(i);
        curvec.save( *(p.cvfileGzStreams[i]), raw_binary );
    }
}