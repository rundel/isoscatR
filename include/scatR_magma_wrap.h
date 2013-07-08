#ifndef _SCATR_MAGMA_WRAP_H
#define _SCATR_MAGMA_WRAP_H

int magma_cholesky( char uplo, int n, double *dA, int ldda, int *info);
void scatr_magma_init();
void scatr_magma_finalize();

#endif
