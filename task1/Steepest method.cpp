#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

//Boost uBLAS
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/lu.hpp"

using namespace boost::numeric::ublas;
int read_answer(char* filename, int N, vector<double>& v)
{
    std::ifstream fin(filename);
    
    if ( !( fin.is_open() ) )
    {
        std::cout << "Failed to open the file " << filename << std::endl;
        return -1; // Failure
    }
    
    for( int i = 0; i < N; i++ )
    {
        fin >> v(i);
    }
    
    return 0;
}
int read_data_from_file(char* filename, int* N, double* E, vector<double>& f, matrix<double>& A )
{
    std::ifstream fin(filename);

    if ( !( fin.is_open() ) )
    {
        std::cout << "Failed to open the file " << filename << std::endl;
        return -1; // Failure
    }

    fin >> (*N);
    fin >> (*E);
    f.resize(*N);
    A.resize(*N, *N);

    for( int i = 0; i < *N; i++ )
    {
        fin >> f(i);
    }

    for( int i = 0; i < *N; i++ )
    {
        for( int j = 0; j < *N; j++ )
        {
            fin >> A(i, j);
        }
    }

    return 0; //success
}

vector<double> calculate( int N, double E, vector<double>& f, matrix<double>& A )
{
    vector<double> u_k (N), u_k_1 (N);
    vector<double> r_k (N);
    vector<double> ans (N);
    //read_answer((char*)"answers.txt", N, ans);
    //std::cout << ans << std::endl;
    double alpha = 0;
    double E_cur = 0;
    int k = 0;

    for( int i = 0; i < N; i++)
    {
        u_k_1(i) = 0;
    }

    do
    {
        u_k = u_k_1;
        r_k = prec_prod(A, u_k) - f;
        
        if (norm_2(r_k) <= E) // Mustn't be a zero value
        {
            break;
        }
        
        alpha = prec_inner_prod(r_k, r_k) / prec_inner_prod(prod(A, r_k), r_k);
        u_k_1 = u_k - alpha * r_k;
        k += 1;
        //E_cur = norm_2(u_k - ans);
        E_cur = norm_2(r_k);

        if (k % 10000 == 0)
	    std::cout  << "E_cur = "<< E_cur << " k = " << k << std::endl;

    } while (E_cur > E);

    std::cout << "Number of iterations: " << k << std::endl;
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
    double E = 0; //condition of termination ||u(k+1) - u(k)|| <= E
    
    vector<double> f(0); //right part of equation
    matrix<double> A(0, 0); //matrix of equation

    if ( read_data_from_file(argv[1], &N, &E, f, A) == -1 )
    {
        exit(EXIT_FAILURE);
    }
    
    f = prec_prod(trans(A), f);
    A = prec_prod(trans(A), A);    

    vector<double> u_final = calculate(N, E, f, A); 
    std::cout << "Epsilon: " << E << std::endl;
    std::cout << "Approximate answer is: " << std::endl;
    
    for (int i = 0; i < N; i++)
        std::cout << std::setprecision(14) << u_final(i) << std::endl;

    return 0;
}
