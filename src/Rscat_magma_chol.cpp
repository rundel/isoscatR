#ifdef USEMAGMA

#include <cublas.h>
#include <cuda.h>
#include <magma.h>
#include "Rscat_magma_chol.h"

void checkCudaError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::stringstream s;
        s << "CUDA " << msg << " : " << cudaGetErrorString(err);
        throw( cuda_exception(s.str()) );
    }
}

void checkCublasError(const char *msg)
{
    cublasStatus err = cublasGetError();
    if(err != CUBLAS_STATUS_SUCCESS) {
        std::stringstream s;
        s << "cuBLAS " << msg << " : " << cublasGetErrorString(err);
        throw( cuda_exception(s.str()) );
    }
}

char *cublasGetErrorString(cublasStatus err)
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


void magChol(arma::mat &B, bool gpu = true)
{
    
   
    int n = a.n_rows, n2 = n * n, info;

    BOOST_ASSET(a.n_rows == a.n_cols);

    double *d_B;
    cublasAlloc(n2, sizeof(double), (void**)&d_B);
    checkCublasError("device memory allocation failed in 'magChol'");
   
    if(gpu) {
      
        cublasSetVector(n2, sizeof(double), B.memptr(), 1, d_B, 1);
        magma_dpotrf_gpu('L', n, d_B, n, &info);
        cublasGetVector(n2, sizeof(double), d_B, 1, B.memptr(), 1);

    } else {
        double *h_B;

        cudaMallocHost((void**)&h_B, n2 * sizeof(double));
        checkCudaError("host memory allocation failed in 'magChol'");

        memcpy(h_B, B, n2 * sizeof(double));
        magma_dpotrf('L', n, h_B, n, &info);
        memcpy(B, h_B, n2 * sizeof(double));

        cudaFreeHost(h_B);
    }
   
    std::stringstream s;
    if (info < 0) {
        s << "illegal argument " << -1 * info << " in 'magChol"; 
        throw( cuda_exception(s.str()) )
    } else if (info > 0) {
        s << "leading minor of order " << info << " is not positive definite";
        throw( cuda_exception(s.str()) )
    }
    cublasFree(d_B);
}

SEXP magma_chol(SEXP rX, SEXP rGPU) {
    
    try {
        Rcpp::NumericMatrix tX(rX);
        arma::mat X = arma::mat(tX.begin(),tX.nrow(),tX.ncol(),false);
        
        bool gpu = Rcpp::as<bool>(rGPU);
        
        magChol(X,gpu);
        
        return Rcpp::wrap(X);
        
    } catch( std::exception &ex ) {
        Rcpp::forward_exception_to_r( ex );
    }
}

#endif