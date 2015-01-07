#ifdef USEMAGMA

#include <iostream>
#include <RcppArmadillo.h>

#include "scatR_gpu.h"
#include "scatR_structs.h"

#include "scatR_magma_wrap.h"

#include <cuda_runtime.h>


void init_GPU_data(GlobalParams &p) {
    
    int n_known = p.nRegions;
    int n_total = p.pred_dist.n_rows;
    int n_pred = n_total-n_known;
    
    scatr_magma_init();

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
    
    scatr_magma_finalize();

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
        Rcpp::Rcout << "Error: CUDA " << msg << " : " << cudaGetErrorString(err) << "\n";
        throw( std::exception() );
    }
}

void checkCublasError(cublasStatus_t err, std::string msg)
{
    std::string err_str;

    switch(err) {
        case CUBLAS_STATUS_SUCCESS :
            err_str = "operation completed successfully";
            break;
        case CUBLAS_STATUS_NOT_INITIALIZED :
            err_str = "CUBLAS library not initialized";
            break;
        case CUBLAS_STATUS_ALLOC_FAILED :
            err_str = "resource allocation failed";
            break;
        case CUBLAS_STATUS_INVALID_VALUE :
            err_str = "unsupported numerical value was passed to function";
            break;
        case CUBLAS_STATUS_ARCH_MISMATCH :
            err_str = "function requires an architectural feature absent from the architecture of the device";
            break;
        case CUBLAS_STATUS_MAPPING_ERROR :
            err_str = "access to GPU memory space failed";
            break;
        case CUBLAS_STATUS_EXECUTION_FAILED :
            err_str = "GPU program failed to execute";
            break;
        case CUBLAS_STATUS_INTERNAL_ERROR :
            err_str = "an internal CUBLAS operation failed";
            break;
        default :
            err_str = "unknown error type";
    }

    if(err != CUBLAS_STATUS_SUCCESS) {
        Rcpp::Rcout << "Error: cuBLAS - " << msg << " - " << err_str;
        throw( std::exception() );
    }
}

#endif