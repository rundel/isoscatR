#ifdef USEMAGMA

#include <iostream>
#include <RcppArmadillo.h>

#include "Rscat_gpu.h"
#include "Rscat_structs.h"

#include <cuda_runtime.h>


GPU_data::GPU_data(GlobalParams const &p, GlobalOptions const &opt) {
    
    int n_known = p.nRegions;
    int n_total = p.pred_dist.n_rows;
    int n_pred = n_total-n_known;
    
    arma::mat dist11 = p.pred_dist(arma::span(0,n_known-1),        arma::span(0,n_known-1));
    arma::mat dist12 = p.pred_dist(arma::span(0,n_known-1),        arma::span(n_known, n_total-1));
    arma::mat dist22 = p.pred_dist(arma::span(n_known, n_total-1), arma::span(n_known, n_total-1));
    
    checkCublasError( cublasCreate_v2(&handle) );
    
    cudaMalloc((void**)&d_dist11, n_known*n_known * sizeof(double));
    cudaMalloc((void**)&d_dist12, n_known*n_pred * sizeof(double));
    cudaMalloc((void**)&d_dist22, n_pred*n_pred * sizeof(double));
    checkCudaError("GPU_data - Init - Dist");
    
    checkCublasError( cublasSetMatrix(n_known, n_known, sizeof(double), dist11.memptr(), n_known, d_dist11, n_known) );
    checkCublasError( cublasSetMatrix(n_known, n_pred,  sizeof(double), dist12.memptr(), n_known, d_dist12, n_known) );
    checkCublasError( cublasSetMatrix(n_pred,  n_pred,  sizeof(double), dist22.memptr(), n_pred,  d_dist22, n_pred) );
    
    cudaMalloc((void**)&d_cov12, n_known*n_pred  * sizeof(double));
    cudaMalloc((void**)&d_cov22, n_pred *n_pred  * sizeof(double));
    cudaMalloc((void**)&d_invcov11, n_known*n_known * sizeof(double));
    cudaMalloc((void**)&d_tmp    , n_pred*n_known * sizeof(double));
    checkCudaError("GPU_data - Init - Misc");
}

GPU_data::~GPU_data(){
    
    cublasDestroy_v2(handle);
    
    cudaFree(d_dist11);
    cudaFree(d_dist12);
    cudaFree(d_dist22);
    
    cudaFree(d_cov12);
    cudaFree(d_cov22);
    
    cudaFree(d_invcov11);
    cudaFree(d_tmp);
    
    checkCudaError("GPU_data - Cleanup");
}


void checkCudaError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::cout << "Error: CUDA " << msg << " : " << cudaGetErrorString(err) << "\n";
        throw( std::exception() );
    }
}

void checkCublasError(cublasStatus_t err)
{
    if(err != CUBLAS_STATUS_SUCCESS) {
        std::cout << "Error: cuBLAS - " << cublasGetErrorString(err);
        throw( std::exception() );
    }
}

std::string cublasGetErrorString(cublasStatus_t err)
{
    switch(err) {
        case CUBLAS_STATUS_SUCCESS :
            return "operation completed successfully";
        case CUBLAS_STATUS_NOT_INITIALIZED :
            return "CUBLAS library not initialized";
        case CUBLAS_STATUS_ALLOC_FAILED :
            return "resource allocation failed";
        case CUBLAS_STATUS_INVALID_VALUE :
            return "unsupported numerical value was passed to function";
        case CUBLAS_STATUS_ARCH_MISMATCH :
            return "function requires an architectural feature absent from the architecture of the device";
        case CUBLAS_STATUS_MAPPING_ERROR :
            return "access to GPU memory space failed";
        case CUBLAS_STATUS_EXECUTION_FAILED :
            return "GPU program failed to execute";
        case CUBLAS_STATUS_INTERNAL_ERROR :
            return "an internal CUBLAS operation failed";
        default :
            return "unknown error type";
    }
}

#endif