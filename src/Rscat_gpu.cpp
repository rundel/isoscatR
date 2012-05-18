#ifdef USEMAGMA

#include <iostream>
#include <RcppArmadillo.h>

#include "Rscat_gpu.h"
#include "Rscat_structs.h"

#include <cuda_runtime.h>


void init_GPU_data(GlobalParams &p) {
    
    int n_known = p.nRegions;
    int n_total = p.pred_dist.n_rows;
    int n_pred = n_total-n_known;
    
    arma::mat dist12 = p.pred_dist(arma::span(0,n_known-1),        arma::span(n_known, n_total-1));
    arma::mat dist22 = p.pred_dist(arma::span(n_known, n_total-1), arma::span(n_known, n_total-1));
    
    checkCublasError( cublasCreate_v2(&p.handle), "handle (Create)" );
    
    cudaMalloc((void**) &p.d_dist12, n_known*n_pred*sizeof(double));  checkCudaError("dist12 (Malloc)");
    checkCublasError( 
        cublasSetMatrix(n_known, n_pred, sizeof(double), dist12.memptr(), n_known, p.d_dist12, n_known), 
        "dist12 (Set)" 
    );
    
    cudaMalloc((void**) &p.d_dist22, n_pred*n_pred * sizeof(double));   checkCudaError("dist22 (Malloc)");
    checkCublasError( 
        cublasSetMatrix(n_pred, n_pred, sizeof(double), dist22.memptr(), n_pred, p.d_dist22, n_pred), 
        "dist22 (Set)"
    );
    
    cudaMalloc((void**) &p.d_cov12,    n_known * n_pred  * sizeof(double)); checkCudaError("cov12 (Malloc)");
    cudaMalloc((void**) &p.d_cov22,    n_pred  * n_pred  * sizeof(double)); checkCudaError("cov22 (Malloc)");
    cudaMalloc((void**) &p.d_invcov11, n_known * n_known * sizeof(double)); checkCudaError("invcov11 (Malloc)");
    cudaMalloc((void**) &p.d_tmp,      n_pred  * n_known * sizeof(double)); checkCudaError("tmp (Malloc)");
}

void cleanup_GPU_data(GlobalParams &p) {
    
    checkCublasError( 
        cublasDestroy_v2(p.handle), "handle (Destroy)" 
    );
    
    cudaFree(p.d_dist12);   checkCudaError("dist12 (Free)");
    cudaFree(p.d_dist22);   checkCudaError("dist22 (Free)");
    cudaFree(p.d_cov12);    checkCudaError("cov12 (Free)");
    cudaFree(p.d_cov22);    checkCudaError("cov22 (Free)");
    cudaFree(p.d_invcov11); checkCudaError("invcov11 (Free)");
    cudaFree(p.d_tmp);      checkCudaError("tmp (Free)");
}


void checkCudaError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::cout << "Error: CUDA " << msg << " : " << cudaGetErrorString(err) << "\n";
        throw( std::exception() );
    }
}

void checkCublasError(cublasStatus_t err, std::string msg)
{
    if(err != CUBLAS_STATUS_SUCCESS) {
        std::cout << "Error: cuBLAS - " << msg << " - " << cublasGetErrorString(err);
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