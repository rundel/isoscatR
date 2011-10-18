#include <cuda.h>

__global__ void powered_exponential_kernel(double* dist, int n, double sigma2, double phi, double kappa, double nugget) {

    const int xid = blockIdx.x * blockDim.x + threadIdx.x,
              yid = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (xid < n && yid < n) {
        dist[yid+n*xid] = sigma2 * exp( -pow(dist[yid+n*xid] / phi, kappa) ) + nugget*(xid==yid);
    }
}



void cov_powered_exponential_gpu(double* dist, int n, double sigma2, double phi, double kappa, double nugget) {

    int BLK_DX = 16, BLK_DY = 16;

    dim3 grids( (n+BLK_DX-1)/BLK_DX, (n+BLK_DY-1)/BLK_DY );
    dim3 threads( BLK_DX, BLK_DY );

    powered_exponential_kernel<<<grids, threads>>>(dist, n, sigma2, phi, kappa, nugget);
    
    cudaThreadSynchronize();
}






