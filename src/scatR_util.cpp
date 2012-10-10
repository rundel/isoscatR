#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "scatR_util.h"

void open_cvfiles(GlobalParams &p, GlobalOptions &opt) {
    
    p.cvfileStreams.resize(p.locate_indivs.size());
    p.cvfileGzStreams.resize(p.locate_indivs.size());

    for(int l=0; l<p.locate_indivs.size(); l++) {
        
        std::string file = opt.TMPDIR + "/CV_Ind" 
                      + boost::lexical_cast<std::string>(p.locate_indivs[l]) + "_"
                      + boost::lexical_cast<std::string>(p.chain_num) + ".gz";
        
        p.cvfileStreams[l] = new std::ofstream(file.c_str(), std::ios_base::binary);
        
        p.cvfileGzStreams[l] = new boost::iostreams::filtering_ostream;
        p.cvfileGzStreams[l]->push( boost::iostreams::gzip_compressor() );
        p.cvfileGzStreams[l]->push( *(p.cvfileStreams[l]) );
    }
}

void close_cvfiles(GlobalParams &p, GlobalOptions &opt) {
    
    for(int l=0; l<p.locate_indivs.size(); l++) {
        delete p.cvfileGzStreams[l];  
        delete p.cvfileStreams[l];
    }
}


void open_allelefiles(GlobalParams &p, GlobalOptions &opt) {
    
    p.alfileStreams.resize(p.nLoci);
    p.alfileGzStreams.resize(p.nLoci);
 
    for(int l=0; l<p.nLoci; l++) {
        
        p.alfileStreams[l].resize(p.nAlleles[l]);
        p.alfileGzStreams[l].resize(p.nAlleles[l]);
        
        for(int j=0; j<p.nAlleles[l]; j++) {
            std::string file = opt.TMPDIR + "/Al" 
                                + boost::lexical_cast<std::string>(l+1) + "-"
                                + boost::lexical_cast<std::string>(j+1) + "_"
                                + boost::lexical_cast<std::string>(p.chain_num) + ".gz";
        
            p.alfileStreams[l][j] = new std::ofstream(file.c_str(), std::ios_base::binary);
        
            p.alfileGzStreams[l][j] = new boost::iostreams::filtering_ostream;
            p.alfileGzStreams[l][j]->push( boost::iostreams::gzip_compressor() );
            p.alfileGzStreams[l][j]->push( *(p.alfileStreams[l][j]) );
        }
    }
}
 
void close_allelefiles(GlobalParams &p, GlobalOptions &opt) {
    
    for(int l=0; l<p.nLoci; l++) {
        for(int j=0; j<p.nAlleles[l]; j++) {
            delete p.alfileGzStreams[l][j];  
            delete p.alfileStreams[l][j];
        }
    }
}

arma::mat calc_rotation_mat(double angle) {
    
    arma::mat res(2,2);
    res(0,0) =  cos(angle);
    res(0,1) = -sin(angle);
    res(1,0) =  sin(angle);
    res(1,1) =  cos(angle);
    
    return res;
}

arma::mat calc_stretch_mat(double ratio) {
    
    arma::mat res(2,2);
    res.eye();
    res(1,1) = 1/ratio;
    
    return res;
}

double calc_multinomial_loglik(const arma::mat &theta, const arma::mat &count, const arma::colvec &sumCount) 
{ 
    return arma::accu( arma::sum(count % theta,1) - sumCount % arma::log(arma::sum(arma::exp(theta),1)) );
}

double calc_multivar_normal_loglik(const arma::mat &x, const arma::mat &mu, const arma::mat &sigmaInv, double sigmaDet) {

    double res = 0.0;
    for(int j=0; j<x.n_cols; j++) {
        arma::mat tmp = -0.5*( log(sigmaDet) + arma::trans(x.col(j) - mu.col(j)) * sigmaInv * (x.col(j) - mu.col(j)) );
        //tmp.print();
        res += arma::accu(tmp);
    }
        
    return res;
}

arma::mat calc_f(const arma::mat &theta) {
    
    int nAlleles = theta.n_cols;
    
    arma::mat expTheta = arma::exp(theta);
    arma::mat sumExpTheta = arma::sum(expTheta,1) * arma::ones<arma::rowvec>(nAlleles);
    
    return expTheta / sumExpTheta;
}


void outputAccepts(GlobalParams &p, GlobalOptions &opt) {
    
    std::cout << std::endl;
    std::cout << "Acceptance Rates:" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    
    std::cout << "alpha  : ";
    for (int l=0; l<p.alpha.size(); l++) 
        std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.alphaAccept[l], p.alphaAttempt[l]) << " ";
    std::cout << std::endl;
    
    if (!opt.FIXANGLE) {
        std::cout << "angle  : ";
        std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.angleAccept, p.angleAttempt) << " ";
        std::cout << std::endl;
    }
    
    if (!opt.FIXRATIO) {
        std::cout << "ratio  : ";
        std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.ratioAccept, p.ratioAttempt) << " ";
        std::cout << std::endl;
    }
    
    if (!opt.FIXXI) {
        std::cout << "xi     : ";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.xiAccept[l],p.xiAttempt[l]) << " ";
        std::cout << std::endl;
    }
    
    if (!opt.FIXBETA) {
        std::cout << "beta   : ";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.betaAccept[l], p.betaAttempt[l]) << " ";
        std::cout << std::endl;
    }
    
    if (!opt.FIXETA) {
        std::cout << "eta    : ";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.etaAccept[l], p.etaAttempt[l]) << " ";
        std::cout << std::endl;
    }
    
    std::cout << "theta  : ";
    for (int l=0; l<p.nLoci; l++) {
        std::cout << std::setprecision(3) << std::setw(5) << accept_ratio(p.thetaAccept[l], p.thetaAttempt[l]) << " ";
    }
    std::cout << std::endl << std::endl;
    
}


void outputTuning(GlobalParams &p, GlobalOptions &opt) {
    
    std::cout << std::endl;
    std::cout << "Tuning Results:" << std::endl;
    std::cout << "=============================================" << std::endl;

    
    for (int l=0; l<opt.ALPHASD.size(); l++) {
        std::cout << "alpha" << l << " : "; 
        std::cout << "[" << std::setprecision(3) << std::setw(5) << opt.ALPHASD(l) << "]";
        std::cout << " " << std::setprecision(3) << std::setw(5) << p.alpha_sd(l) << std::endl;
    }
    
    if (!opt.FIXANGLE) {
        std::cout << "angle  : ";
        std::cout << "[" << std::setprecision(3) << std::setw(5) << opt.ANGLESD << "]";
        std::cout << " " << std::setprecision(3) << std::setw(5) << p.angle_sd << std::endl;
    }
    
    if (!opt.FIXRATIO) {
        std::cout << "ratio  : ";
        std::cout << "[" << std::setprecision(3) << std::setw(5) << opt.RATIOSD << "]";
        std::cout << " " << std::setprecision(3) << std::setw(5) << p.ratio_sd << std::endl;
    }
    
    if (!opt.FIXXI) {
        std::cout << "xi     : [" << std::setprecision(3) << std::setw(5) << opt.XISD << "] ";
        std::cout << "(" << std::setprecision(3) << std::setw(5) << mean(p.xi_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << " " << std::setprecision(3) << std::setw(5) << p.xi_sd(l);
        std::cout << std::endl;
    }
    
    if (!opt.FIXBETA) {
        std::cout << "beta   : [" << std::setprecision(3) << std::setw(5) << opt.BETASD << "] ";
        std::cout << "(" << std::setprecision(3) << std::setw(5) << mean(p.beta_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << " " << std::setprecision(3) << std::setw(5) << p.beta_sd(l);
        std::cout << std::endl;
    }
    
    if (!opt.FIXETA) {
        std::cout << "eta    : [" << std::setprecision(3) << std::setw(5) << opt.ETASD << "] ";
        std::cout << "(" << std::setprecision(3) << std::setw(5) << mean(p.eta_sd) << ")";
        for (int l=0; l<p.nLoci; l++) 
            std::cout << " " << std::setprecision(3) << std::setw(5) << p.eta_sd(l);
        std::cout << std::endl;
    }
    
    std::cout << "theta  : [" << std::setprecision(3) << std::setw(5) << opt.THETASD << "] ";
    std::cout << "(" << std::setprecision(3) << std::setw(5) << mean(p.theta_sd) << ")";
    for(int l=0; l<p.nLoci; l++) {
        std::cout << " " << std::setprecision(3) << std::setw(5) << p.theta_sd[l];
    }
    std::cout << std::endl << std::endl;
}



// compute great circle distance from positions in radians
// from http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/TSPFAQ.html

double great_circle_dist(double x1, double y1, double x2, double y2)
{
    double RRR = 6378.388; 
    
    double q1 = cos( x1 - x2 ); 
    double q2 = cos( y1 - y2 ); 
    double q3 = cos( y1 + y2 ); 
    
    return( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) );
}

double accept_ratio(unsigned int accept, unsigned int attempt)
{
    return (1.0*accept) / attempt;
}


arma::mat distance_mat(arma::mat const &m) 
{
    return distance_mat(m.col(0),m.col(1));
}

arma::mat distance_mat(arma::colvec const &xc, arma::colvec const &yc) {

    arma::mat dist = arma::zeros<arma::mat>(xc.n_elem, xc.n_elem);
    for(int i=1; i < xc.n_elem; ++i) {
        for(int j=0; j < i; ++j) {
            dist(i,j) = great_circle_dist(xc(i), yc(i), xc(j), yc(j));
            dist(j,i) = dist(i,j);
        }
    }
    
    return(dist);
}

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2)
{
    return( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
}

//  Winding number test for a point in a polygon
//  modelled on code from softsurfer.com, by Dan Sunday
//      Input:   x,y = a point,
//               BoundaryX and BoundaryY = points of a polygon with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if (x,y) is outside polygon)

int isInsideBoundary( double x, double y, arma::mat boundary)
{
    if(boundary.n_elem == 0) // if no boundary, just return 1
        return 1;
    
    int wn = 0;    // the winding number counter
    
    for (int i=0; i < (boundary.n_rows-1); i++) {   // edge from V[i] to V[i+1]
        if (boundary(i,1) <= y) {                   // start y <= P.y
            if (boundary(i+1,1) > y 
                && isLeft( boundary(i,0),boundary(i,1), boundary(i+1,0),boundary(i+1,1), x, y) > 0)  
                ++wn;            
        } else {       
            if (boundary(i+1,1) <= y 
                && isLeft( boundary(i,0),boundary(i,1), boundary(i+1,0),boundary(i+1,1), x, y) < 0)  
                --wn;            
        }
    }
    return wn;
}