#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    Define left part of your equation here:
    y''(x) = left_aprt(x, y, y')
*/
double left_part(double y, double dy_dx, double x)
{
    return -dy_dx + y / x + (x + 1) / x;
}

/*
    Returns two arrays of y(x) and y'(x)
    according to the chosen partition of
    segment - [x0; x_right]
    
    alpha - an arbitrarily chosen value of y'(x0)
    N - parameter of partition
    h - step of partition
*/
double** Euler_solve(int N, double h, double x0, double alpha)
{
    double* y1 = (double*) malloc (sizeof(double) * (N + 1));
    double* y2 = (double*) malloc (sizeof(double) * (N + 1));
    double** y = (double**) malloc (sizeof(double*) * 2);

    // Initial conditions
    y[0] = y1;
    y[1] = y2;
    y1[0] = 0;
    y2[0] = alpha;

    // Calculating y and y' in every point of the net(partition)
    for(int i = 1; i <= N; i++)
    {
        y1[i] = y1[i-1] + y2[i-1] * h;
        y2[i] = y2[i-1] + left_part(y1[i-1], y2[i-1], x0 + (i - 1) * h ) * h;
    }

    return y;
}
