#include "Rscat.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

SEXP read_allele_file(SEXP Rfile) {
    
    string file_name = as<string>(Rfile);
    
    ifstream filestream(file_name.c_str(), ios_base::in | ios_base::binary);
    
    
    if (filestream.is_open()) {
        
        boost::iostreams::filtering_istream f;
        f.push( boost::iostreams::gzip_decompressor() );
        f.push( filestream );
        
        vector<char> data;
        boost::iostreams::copy(f, boost::iostreams::back_inserter(data));
        
        vec x( (double *) &data[0], data.size() / sizeof(double), false);
        
        return(wrap(x));

    } else {
        cout << "Unable to open allele file (" << file_name <<")" << endl;
        return(R_NilValue);
    }
}


SEXP prec_sum(SEXP Rvec) {

    // Full precision summation using multiple floats for intermediate values
    // Rounded x+y stored in hi with the round-off stored in lo.  Together
    // hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
    // to each partial so that the list of partial sums remains exact.
    // Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
    // www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps

    NumericVector vec = as<NumericVector>(Rvec);

    if (vec.size() == 1)
        return(Rvec);

    vector<double> partials(1,vec[0]);               // sorted, non-overlapping partial sums
    int k;

    for (int i=1; i<vec.size(); i++) {
        double x = vec[i];
        k = 0;

        for (int j=0; j<partials.size(); j++) {

            double y = partials[j];

            if (abs(x) < abs(y)) {
                double tmp=y;
                y=x;
                x=tmp;
            }

            double hi = x + y;
            double lo = y - (hi - x);
            if (lo != 0.0) {
                partials[k] = lo;
                k++;
            }
            x = hi;
        }

        partials.resize(k+1);
        partials[k] = x;
    }

    double sum = 0.0;
    for (int j=0; j<partials.size(); j++) 
        sum += partials[j];

    return (wrap(sum));
}

mat calc_rotation_mat(double angle) {
    
    mat res(2,2);
    res(0,0) =  cos(angle);
    res(0,1) = -sin(angle);
    res(1,0) =  sin(angle);
    res(1,1) =  cos(angle);
    
    return(res);
}

mat calc_stretch_mat(double ratio) {
    
    mat res(2,2);
    res.eye();
    res(1,1) = 1/ratio;
    
    return(res);
}

colvec calc_multinom_loglik( const mat theta, const mat count, const colvec sumCount) {
    
    //mat expTheta = exp(theta);
    //mat invLogit = expTheta / (1+expTheta);
    //mat sumInvLogit = sum(invLogit,1);
    //
    //return( sum(count % log(invLogit),1) - sumCount % log(sumInvLogit) );

    return( sum(count % theta,1) - sumCount % log(sum(exp(theta),1)) );
}

mat calc_theta(const double &mu, const rowvec &eta, const double &xi, const mat &L, const mat &X) {
    int nRegions = X.n_rows;

    return( ones<colvec>(nRegions) * (mu + xi * eta) + L * X );
}


mat calc_f(const mat &theta) {
    
    int nAlleles = theta.n_cols;
    
    mat expTheta = exp(theta);
    mat sumExpTheta = sum(expTheta,1) * ones<rowvec>(nAlleles);
    
    return( expTheta / sumExpTheta );
}


void parseOptions(SEXP sexpRopt, GlobalOptions &opt) {
    List Ropt(sexpRopt);
    
    opt.VERBOSE        = as<bool>(Ropt["VERBOSE"]);
    
    opt.TMPDIR         = as<string>(Ropt["TMPDIR"]);
    opt.FILESUFFIX     = as<string>(Ropt["FILESUFFIX"]);
    
    opt.ADAPT          = as<bool>(Ropt["ADAPT"]);
    opt.TUNEINTERVAL    = as<int>(Ropt["TUNEINTERVAL"]);

    opt.LOCATE         = as<bool>(Ropt["LOCATE"]);
    opt.MAXCELL        = as<double>(Ropt["MAXCELL"]);
    opt.CROSSVALIDATE  = as<bool>(Ropt["CROSSVALIDATE"]);
    
    
    opt.RETURNFIT    = as<bool>(Ropt["RETURNFIT"]);
    opt.USEMATERN      = as<bool>(Ropt["USEMATERN"]);
 
    opt.PSEUDOCOUNT    = as<double>(Ropt["PSEUDOCOUNT"]);

    opt.FIXALPHA       = as<LogicalVector>(Ropt["FIXALPHA"]);
    opt.ALPHA          = as<NumericVector>(Ropt["ALPHA"]);
    opt.ALPHAMIN       = as<NumericVector>(Ropt["ALPHAMIN"]);
    opt.ALPHAMAX       = as<NumericVector>(Ropt["ALPHAMAX"]);

    opt.ALPHASD        = as<NumericVector>(Ropt["ALPHASD"]);
    
    opt.ANGLE          = as<double>(Ropt["ANGLE"]);
    opt.FIXANGLE       = as<bool>(Ropt["FIXANGLE"]);
    opt.ANGLESD        = as<double>(Ropt["ANGLESD"]);
    
    opt.RATIO          = as<double>(Ropt["RATIO"]);
    opt.FIXRATIO       = as<bool>(Ropt["FIXRATIO"]);
    opt.RATIOSD        = as<double>(Ropt["RATIOSD"]);
    
    opt.XIRANGE        = as<NumericVector>(Ropt["XIRANGE"]);
    opt.FIXXI          = as<bool>(Ropt["FIXXI"]);
    opt.XI             = as<NumericVector>(Ropt["XI"]);
    opt.XISD           = as<double>(Ropt["XISD"]);
    opt.SIGMAXI        = as<double>(Ropt["SIGMAXI"]);
    
    opt.FIXMU          = as<bool>(Ropt["FIXMU"]);
    opt.MU             = as<NumericVector>(Ropt["MU"]);
    opt.MUSD           = as<double>(Ropt["MUSD"]);
    
    opt.FIXETA         = as<bool>(Ropt["FIXETA"]);
    opt.ETA            = as<NumericVector>(Ropt["ETA"]);
    opt.ETASD          = as<double>(Ropt["ETASD"]);
    
    opt.BETARANGE      = as<NumericVector>(Ropt["BETARANGE"]);
    opt.FIXBETA        = as<bool>(Ropt["FIXBETA"]);
    opt.BETA           = as<NumericVector>(Ropt["BETA"]);
    opt.BETASD         = as<double>(Ropt["BETASD"]);
    opt.SIGMABETA      = as<double>(Ropt["SIGMABETA"]);
    
    
    opt.XSD            = as<double>(Ropt["XSD"]);
    
    opt.LOCALSD        = as<double>(Ropt["LOCALSD"]);
    opt.GLOBALSD       = as<double>(Ropt["GLOBALSD"]);
    
    opt.NULLPROB       = as<double>(Ropt["NULLPROB"]);
    opt.DELTA          = as<double>(Ropt["DELTA"]);
    
    return;
}



// compute great circle distance from positions in radians
//from http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/TSPFAQ.html
double calc_distance(double x1, double y1, double x2, double y2) {
    double RRR = 6378.388; 
    
    double q1 = cos( x1 - x2 ); 
    double q2 = cos( y1 - y2 ); 
    double q3 = cos( y1 + y2 ); 
    
    return( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) );
}

SEXP R_calc_distance_to_point(SEXP px, SEXP py, SEXP x, SEXP y) {
    NumericVector xc(x), yc(y);
    
    NumericVector dist(xc.size());
    for(int i=0; i < xc.size(); ++i) {
        dist[i] = calc_distance(as<double>(px), as<double>(py), xc[i], yc[i]);
    }
    
    return(dist);
}


double calc_accept_ratio(unsigned int accept, unsigned int attempt) {
    
    return( (1.0*accept) / attempt );
}


mat calc_distance_mat(mat const &m) {
    return(calc_distance_mat(m.col(0),m.col(1)));
}

mat calc_distance_mat(colvec const &xc, colvec const &yc) {

    mat dist = zeros<mat>(xc.n_elem, xc.n_elem);
    for(int i=1; i < xc.n_elem; ++i) {
        for(int j=0; j < i; ++j) {
            dist(i,j) = calc_distance(xc(i), yc(i), xc(j), yc(j));
            dist(j,i) = dist(i,j);
        }
    }
    
    return(dist);
}

SEXP R_calc_distance(SEXP x, SEXP y) {
    NumericVector xt(x);                 // creates Rcpp vector from SEXP
    NumericVector yt(y);                   // creates Rcpp matrix from SEXP
    
    colvec xc(xt.begin(), xt.size(), false);
    colvec yc(yt.begin(), yt.size(), false);
    
    mat res = calc_distance_mat(xc,yc);
    
    return(wrap(res));
}

// normal density
colvec dnorm(colvec x) {
    return ((1.0/sqrt(2*math::pi())) * exp((-0.5)*square(x)));
}
double dnorm(double x) {
    return ((1.0/sqrt(2*math::pi())) * exp((-0.5)*(x*x)));
}


double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2) {
    return( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
}

//  Winding number test for a point in a polygon
//  modelled on code from softsurfer.com, by Dan Sunday
//      Input:   x,y = a point,
//               BoundaryX and BoundaryY = points of a polygon with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if (x,y) is outside polygon)
int isInsideBoundary( double x, double y, mat boundary) {
    
    if(boundary.n_elem == 0) // if no boundary, just return 1
        return 1;
  

    int wn = 0;    // the winding number counter
    
    // loop through all edges of the polygon
    for (int i=0; i < (boundary.n_rows-1); i++) {   // edge from V[i] to V[i+1]
        if (boundary(i,1) <= y) {         // start y <= P.y
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