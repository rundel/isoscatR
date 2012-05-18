#ifndef _RSCAT_GPU_H
#define _RSCAT_GPU_H

#include <cublas_api.h>
#include "Rscat_structs.h"

void init_GPU_data(GlobalParams &p);
void cleanup_GPU_data(GlobalParams &p);

void cov_powered_exponential_gpu(double* dist, double* cov,
                                 const int n, const int m,
                                 double sigma2, double phi, 
                                 double kappa, double nugget,
                                 int n_threads);


void checkCudaError(const char *msg);
void checkCublasError(cublasStatus_t err, std::string msg);
std::string cublasGetErrorString(cublasStatus_t err);



#endif