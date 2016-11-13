#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int read_data_from_file(char* filename, int* N, long double* E, long double** f, long double*** A )
{
    std::ifstream fin(filename);

    if ( !( fin.is_open() ) )
    {
        std::cout << "Failed to open the file " << filename << std::endl;
        return -1; // Failure
    }
        
    fin >> (*N);
    fin >> (*E);
    
    (*f) = (long double*) malloc ( (*N) * sizeof(long double) );
    (*A) = (long double**) malloc ( (*N) * sizeof(long double*) );
    
    if ( f == NULL || A == NULL )
    {
        std::cout << "Null pointer was given to the function" << std::endl;
        return -1;
    }

    for( int i = 0; i < (*N); i++ )
    {
        fin >> (*f)[i];
    }
    
    for( int i = 0; i < (*N); i++ )
    {
        (*A)[i] = (long double*) malloc ( (*N) * sizeof(long double) );
        
        if ( (*A)[i] == NULL )
        {
            std::cout << "Null pointer was given to the function" << std::endl;
            return -1;
        }
        
        for( int j = 0; j < (*N); j++ )
        {
            
            fin >> (*A)[i][j];
        }
    }

    return 0; //success
}

long double* M_v_mul( int N, long double** M, long double* v)
{
    long double* ans = (long double*) malloc ( N * sizeof(long double) );
    
    for( int i = 0; i < N; i++ )
    {
        ans[i] = 0;
        
        for( int j = 0; j < N; j++)
        {
            ans[i] += M[i][j] * v[j];
        }
    }
    
    return ans;
}

long double scalar_mul( int N, long double* v1, long double* v2)
{
    long double ans = 0;
    
    for( int i = 0; i < N; i++ )
    {   
        ans += v1[i] * v2[i];
    }
    
    return ans;
}

long double* s_v_mul( int N, long double s, long double* v, long double* d)
{
    for( int i = 0; i < N; i++ )
    {   
        d[i] = s * v[i];
    }
    
    return d;
}

long double* sum( int N, long double* v1, long double* v2, long double* ans, int coeff)
{   
    for( int i = 0; i < N; i++ )
    {   
        ans[i] = v1[i] + coeff * v2[i];
    }
    
    return ans;
}

long double norm(int N, long double* v )
{
    long double sum = 0;
    
    for ( int i = 0; i < N; i++)
    {
        sum += abs(v[i]);// * v[i];
    }
    
    long double ans = sum;//sqrtl( sum );
    return ans;
}

long double* calculate( int N, long double E, long double* f, long double** A )
{
    long double* u_k = (long double*) calloc ( N, sizeof(long double) );
    long double* u_k_1 = (long double*) calloc ( N, sizeof(long double) );
    long double* tmp_alloc = (long double*) malloc ( N * sizeof(long double) );
    long double* r_k = NULL;
    long double* tmp = NULL;
    long double alpha = 0;
    int k = 0;
    long double E_cur = 0;

    do
    {
        tmp = u_k;
        u_k = u_k_1;
        u_k_1 = tmp;
        
        free(r_k);
        r_k = NULL;
        
        r_k = M_v_mul(N, A, u_k);
        
        sum( N, r_k, f, r_k, -1);
        
        if (norm( N, r_k ) <= E)
        {            
            tmp = u_k;
            u_k = u_k_1;
            u_k_1 = tmp;
            
            break;
        }
        
        tmp = M_v_mul( N, A, r_k );
        alpha = scalar_mul( N, r_k, r_k ) / scalar_mul( N, tmp, r_k );
        
        free(tmp);
        tmp = NULL;
        
        s_v_mul( N, alpha, r_k, tmp_alloc );
        sum( N, u_k, tmp_alloc, u_k_1, -1 );
        
        k += 1;
        
        E_cur = norm( N, r_k );
        
        if ( k % 10000 == 0 )
        {
            std::cout << "E_cur: " << E_cur << std::endl << "K: "<< k << std::endl;
        }
    } while (E_cur > E);

    std::cout << "Number of iterations: " << k << std::endl;
    
    free(u_k);
    u_k = NULL;
    
    free(tmp_alloc);
    tmp_alloc = NULL;
    
    free(r_k);
    r_k = NULL;
    
    return u_k_1;
}

int main( int argc, char** argv )
{
    if ( argc != 2)
    {
        std::cout << "You should enter THE ONLY parameter == filename" << std::endl;
        exit(EXIT_FAILURE);
    }
 
    int N = 0; //dimension
    long double E = 0; //condition of termination ||u(k+1) - u(k)|| <= E
    long double* f  = NULL;// right part of equation
    long double** A = NULL; // matrix of equation

    if ( read_data_from_file(argv[1], &N, &E, &f, &A) == -1 )
    {
        exit(EXIT_FAILURE);
    }
    
    long double* u_final = calculate(N, E, f, A);
    std::cout << "Epsilon: " << E << std::endl;
    std::cout << "Approximate answer is: " << std::endl;
    
    for (int i = 0; i < N; i++)
        std::cout << std::setprecision(5) << u_final[i] << std::endl;
/*    
    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            std::cout << A[i][j] <<" ";
        }
        std::cout << std::endl;
    }
*/
    free(f);
    free(u_final);
    
    for(int i = 0; i < N; i++)
    {
        free(A[i]);
    }
    free(A);
    
    return 0;
}