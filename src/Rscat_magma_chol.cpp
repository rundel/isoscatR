#ifdef USEMAGMA

#include "Rscat.h"
#include <magma.h>
#include <boost/timer.hpp>

SEXP magma_chol(SEXP rX, SEXP rGPU, SEXP rFLOAT)
{    
    arma::mat X = Rcpp::as<arma::mat>(rX);

    bool gpu = Rcpp::as<bool>(rGPU);
    
    mag_schol_debug(X,gpu);
    
    return Rcpp::wrap(X);
}

SEXP pow_exp_test(SEXP rX, SEXP ralpha, SEXP gpu, SEXP threads)
{    
    arma::mat dist = Rcpp::as<arma::mat>(rX);
    Rcpp::NumericVector alpha(ralpha); 
    
    if (Rcpp::as<bool>(gpu)) {
        int n = dist.n_rows, n2=n*n;
        double *d_B;
    
        cublasAlloc(n2, sizeof(double), (void**)&d_B);
        cublasSetVector(n2, sizeof(double), dist.memptr(), 1, d_B, 1);        

        cov_powered_exponential_gpu(d_B,n,alpha[0],alpha[1],alpha[2],alpha[3],Rcpp::as<int>(threads));
        checkCudaError("Cov Power Exp Kernel");
    
        cublasGetVector(n2, sizeof(double), d_B, 1, dist.memptr(), 1);
        cublasFree(d_B);
        
        return Rcpp::wrap(dist);
    } else {
        return Rcpp::wrap( cov_powered_exponential(alpha[0], alpha[1], alpha[2], alpha[3], dist) );
    }
}




arma::mat calc_L_gpu( std::vector<double> alpha, arma::mat& dist, bool usematern) {
    
    if (usematern)
        return calc_L(alpha, dist, usematern);
    
    int n = dist.n_rows, n2 = n*n, info;
    double* d_B;
    
    cublasAlloc(n2, sizeof(double), (void**)&d_B);
    cublasSetVector(n2, sizeof(double), dist.memptr(), 1, d_B, 1);        
    
    cov_powered_exponential_gpu(d_B,n,alpha[0],alpha[1],alpha[2],alpha[3],64);
    checkCudaError("Cov Power Exp Kernel");
    
    magma_dpotrf_gpu('L', n, d_B, n, &info);
   
    if (info < 0) {
        std::cout << "Error: illegal argument " << -1 * info << " in magChol\n"; 
        throw( std::exception() );
    } else if (info > 0) {
        std::cout << "Error: leading minor of order " << info << " is not positive definite\n";
        throw( std::exception() );
    }
    
    arma::mat L(dist.n_rows, dist.n_cols);
    cublasGetVector(n2, sizeof(double), d_B, 1, L.memptr(), 1);
    cublasFree(d_B);
    
    for(int col=1; col < n; ++col) {
        double* colptr = L.colptr(col);
        memset(colptr, 0, col * sizeof(double));
    }
    
    return(L);
}


void mag_schol_debug(arma::mat &A, bool gpu)
{
    boost::timer t,total;
    arma::fmat B = arma::conv_to<arma::fmat>::from(A);
    std::cout << "Conv: " << t.elapsed() << "\n"; 
    
    int n = B.n_rows, n2 = n * n, info;
    
    BOOST_ASSERT(B.n_rows == B.n_cols);
    
    if(gpu) {
        float *d_B;
        t.restart();
        cublasAlloc(n2, sizeof(float), (void**)&d_B);
        checkCublasError("device memory allocation failed in 'magChol'");
        std::cout << "Alloc: " << t.elapsed() << "\n"; 
        
        t.restart();
        cublasSetVector(n2, sizeof(float), B.memptr(), 1, d_B, 1);
        std::cout << "Set: " << t.elapsed() << "\n"; 
        
        t.restart();
        magma_spotrf_gpu('L', n, d_B, n, &info);
        std::cout << "Calc: " << t.elapsed() << "\n"; 
        
        t.restart();
        cublasGetVector(n2, sizeof(float), d_B, 1, B.memptr(), 1);
        std::cout << "Get: " << t.elapsed() << "\n"; 
        
        t.restart();
        cublasFree(d_B);
        std::cout << "Free: " << t.elapsed() << "\n"; 
        
    } else {
        float *h_B;
        cudaMallocHost((void**)&h_B, n2 * sizeof(float));
        checkCudaError("host memory allocation failed in 'magChol'");

        memcpy(h_B, B.memptr(), n2 * sizeof(float));
        magma_spotrf('L', n, h_B, n, &info);
        memcpy(B.memptr(), h_B, n2 * sizeof(float));

        cudaFreeHost(h_B);
    }
   
    if (info < 0) {
        std::cout << "Error: illegal argument " << -1 * info << " in magChol\n"; 
        throw( std::exception() );
    } else if (info > 0) {
        std::cout << "Error: leading minor of order " << info << " is not positive definite\n";
        throw( std::exception() );
    }
    
    t.restart();    
    for(int col=1; col < n; ++col) {
        float* colptr = B.colptr(col);
        memset(colptr, 0, col * sizeof(float));
    }
    std::cout << "ToLow: " << t.elapsed() << "\n"; 
    
    
    t.restart();    
    A = arma::conv_to<arma::mat>::from(B);
    std::cout << "Conv: " << t.elapsed() << "\n"; 
    
    std::cout << "Total: " << total.elapsed() << "\n\n";     
}

void checkCudaError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        std::cout << "Error: CUDA " << msg << " : " << cudaGetErrorString(err) << "\n";
        throw( std::exception() );
    }
}

void checkCublasError(const char *msg)
{
    cublasStatus err = cublasGetError();
    if(err != CUBLAS_STATUS_SUCCESS) {
        std::stringstream s;
        std::cout << "Error: cuBLAS " << msg << " : " << cublasGetErrorString(err);
        throw( std::exception() );
    }
}

std::string cublasGetErrorString(cublasStatus err)
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