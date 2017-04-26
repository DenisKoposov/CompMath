#include <stdio.h>
#include <stdlib.h>

/***
    N - parameter of splitting(N+1 equations)
    a, b, c, d - sets of (N+1) coefficients
***/
double* TDMA_solve(int N, double* a, double* b, double* c, double* d)
{
    double* alpha = (double*) malloc (sizeof(double) * N);
    double* beta  = (double*) malloc (sizeof(double) * N);
    double* y     = (double*) malloc (sizeof(double) * (N + 1));

    alpha[0] = -a[0]/b[0]; 
    beta[0]  = d[0]/b[0];
    // Forward propagation
    for (int i = 1; i < N; i++)
    {
        alpha[i] = -a[i] / (b[i] + c[i] * alpha[i-1]);
        beta[i]  = (d[i] - c[i] * beta[i-1]) / (b[i] + c[i] * alpha[i-1]);
    }

    y[N] = (d[N] - c[N] * beta[N-1]) / (b[N] + c[N] * alpha[N-1]);
    // Backward propagation
    for (int i = N-1; i >= 1; i--)
    {
        y[i] = alpha[i] * y[i+1] + beta[i];
    }

    free(alpha);
    alpha = NULL;
    free(beta);
    beta = NULL;

    return y;
}
