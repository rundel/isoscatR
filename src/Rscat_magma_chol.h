#ifndef _RSCAT_MAGMA_CHOL_H
#define _RSCAT_MAGMA_CHOL_H

#include <cublas.h>
#include <cuda.h>
#include <magma.h>

struct cuda_exception : public std::exception
{
    std::string s;
    cuda_exception(std::string ss)
      : s(ss) {}
   
    const char* what() const throw() 
    { 
        return s.c_str();
    }
};

void checkCudaError(const char *msg);
void checkCublasError(const char *msg);
char *cublasGetErrorString(cublasStatus err);
void magChol(arma::mat &B, bool gpu = true);

SEXP magma_chol(SEXP rX, SEXP rGPU);

#endif