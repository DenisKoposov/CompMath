#include <stdio.h>
#include <stdlib.h>
#include "TDMA.c"
#include "euler.c"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

double p(double x) { return 1; }
double q(double x) { return -1/x; }
double f(double x) { return (x + 1) / x; }

double** delete_answer(double** ans)
{
    free(ans[0]);
    free(ans[1]);
    free(ans);
    return NULL;
}

int main(int argc, char** argv)
{
    // Parameters of partition(net)
    int N = 1000;
    double x0 = 0.5;
    double x_right = 1.0;
    double h = (x_right - x0) / N;
    
/*
    TDMA - tridiagonal matrix algorithm for the following problem:

    y''(x) + p(x)y'(x) + q(x)y(x) = f(x), x0 <= x <= x_right
    w0 * y(x0) + w1 * y'(x0) = g0
    m0 * y(x_right) + m1 * y'(x_right) = g1
*/
    double* a = (double*) malloc (sizeof(double) * (N+1));
    double* b = (double*) malloc (sizeof(double) * (N+1));
    double* c = (double*) malloc (sizeof(double) * (N+1));
    double* d = (double*) malloc (sizeof(double) * (N+1));

    double w0 = 1.0;
    double w1 = 0;
    double g0 = 0;

    double m0 = 0;
    double m1 = 1.0;
    double g1 = 1 - log(0.5);
    
    a[0] = w1;
    b[0] = w0 * h - w1;
    c[0] = 0;
    d[0] = g0 * h;
    
    // Calculating matrix coefficients with recursive relations
    for (int i = 1; i < N; i++)
    {
        double x = x0 + i * h;
        
        a[i] = 2 + p(x) * h;
        b[i] = -4 + 2 * q(x) * h * h;
        c[i] = 2 - p(x) * h;
        d[i] = 2 * f(x) * h * h;
    }

    a[N] = 0;
    b[N] = m0 * h + m1;
    c[N] = -m1;
    d[N] = g1 * h;
    
    double* y = TDMA_solve(N, a, b, c, d);
    FILE* output = fopen("data_TDMA.csv", "w");
    
    for (int i = 0; i <= N; i++)
    {
        fprintf(output, "%f,%f\n", x0 + i * h, y[i]);
    }

    fclose(output);
    free(a);
    a = NULL;
    free(b);
    b = NULL;
    free(c);
    c = NULL;
    free(d);
    d = NULL;
    free(y);
    y = NULL;

    // Shooting method
    double alpha1 = 0; // y'(alpha1, x_right) > y_initial'(x_right)
    int flag1 = 0;
    double alpha2 = 0; // y'(alpha2, x_right) < y_initial'(x_right)
    int flag2 = 0;
    double eps = 0.001; // possible deviation from y_initial(x_right)
    int i = 0;

    while (!flag1)
    {

        double** ans = Euler_solve(N, h, x0,  i);
        if (ans[1][N] >= g1)
        {
            alpha1 = i;
            flag1 = 1;
        }
        else
        {
            alpha2 = i;
            flag2 = 1;
        }

        i += 1;
        delete_answer(ans);
    }

    i = -1;

    while (!flag2)
    {
        double** ans = Euler_solve(N, h, x0, i);
        if (ans[1][N] <= g1)
        {
            alpha2 = i;
            flag2 = 1;
        }

        i -= 1;
        delete_answer(ans);
    }

    double** ans = Euler_solve(N, h, x0, (alpha1 + alpha2) / 2);
    while (fabs(ans[1][N] - g1) > eps)
    {
        if (ans[1][N] > g1)
            alpha1 = (alpha1 + alpha2) / 2;
        else
            alpha2 = (alpha1 + alpha2) / 2;

        delete_answer(ans);
        ans = Euler_solve(N, h, x0, (alpha1 + alpha2) / 2);
    }

    output = fopen("data_SM.csv", "w");

    for (int i = 0; i <= N; i++)
    {
        fprintf(output, "%f,%f,%f\n", x0 + i * h, ans[0][i], ans[1][i]);
    }

    fclose(output);
    delete_answer(ans);

    return 0;
}
