#ifdef USEMAGMA

#include <cuda_runtime_api.h>
#include <cublas.h>
#include <magma.h>

int magma_cholesky( char uplo, int n, double *dA, int ldda, int *info)
{
    return magma_dpotrf_gpu( uplo, n, dA, ldda, info);
}

void scatr_magma_init()
{
    magma_init();
}

void scatr_magma_finalize()
{
    magma_finalize();
}


#endif
