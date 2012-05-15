__global__ void powered_exponential_kernel(double* dist, double* cov,
                                           const int n, const int nm, 
                                           const double sigma2, const double phi, 
                                           const double kappa, const double nugget) 
{
    int n_threads = gridDim.x * blockDim.x;
    int pos = blockDim.x * blockIdx.x + threadIdx.x;

    for (int i = pos; i < nm; i += n_threads)
        cov[i] = sigma2 * exp( -pow(dist[i] / phi, kappa) ) + nugget*( i%n == 0 );
    
}


void cov_powered_exponential_gpu(double* dist, double* cov,
                                 const int n, const int m,
                                 double sigma2, double phi, 
                                 double kappa, double nugget,
                                 int n_threads) 
{
    int nm = n*m;
    int blocks = (n+n_threads-1)/n_threads;
    
    powered_exponential_kernel<<<blocks, n_threads>>>(dist, cov, n, nm, sigma2, phi, kappa, nugget);
}
