#ifndef _RSCAT_GPU_H
#define _RSCAT_GPU_H

#include <cublas_api.h>

//#include "Rscat_structs.h"

void cov_powered_exponential_gpu(double* dist, double* cov,
                                 const int n, const int m,
                                 double sigma2, double phi, 
                                 double kappa, double nugget,
                                 int n_threads);


void checkCudaError(const char *msg);
void checkCublasError(cublasStatus_t err);
std::string cublasGetErrorString(cublasStatus_t err);

class GlobalParams;
class GlobalOptions;

class GPU_data {
  public:
    double *d_invcov11;
    double *d_dist11, *d_dist12, *d_dist22;
    double *d_cov11,  *d_cov12,  *d_cov22;
    double *d_tmp;
    double *d_predcov, *d_predmean;
    cublasHandle_t handle;
    
    GPU_data() {}
    GPU_data(GlobalParams const &p, GlobalOptions const &opt);
    ~GPU_data();
};

#endif