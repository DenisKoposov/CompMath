#include <iostream>
#include <fstream>
#include <cstdlib>

//Boost uBLAS
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
using namespace boost::numeric::ublas;


void write_data_to_file(const char* filename, int N, double E, matrix<int>& A, vector<int>& f)
{
    std::ofstream fout(filename);
    
    if ( !( fout.is_open() ) )
    {
        std::cout << "Failed to open the file " << filename << std::endl;
        exit(EXIT_FAILURE); // Failure
    }

    fout << N << std::endl;
    fout << E << std::endl;
    
    for( int i = 0; i < N; i++ )
    {
        fout << f(i) << std::endl;
    }

    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            fout << A(i, j) << std::endl;
        }
    }
}

void write_answer_to_file(const char* filename, vector<int>& x)
{
    std::ofstream fout(filename);
    
    if ( !( fout.is_open() ) )
    {
        std::cout << "Failed to open the file " << filename << std::endl;
        exit(EXIT_FAILURE); // Failure
    }

    for( int i = 0; i < x.size(); i++ )
    {
        fout << x(i) << std::endl;
    }
}

int main( int argc, char** argv )
{
    if (argc != 3)
    {
        std::cout << "Incorrect number of parameters" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    int N = strtoll(argv[1], NULL, 10); //dimension
    double E = strtod(argv[2], NULL); //condition of termination ||u(k+1) - u(k)|| <= E
    
    vector<int> x(N), f(N); //right part of equation
    matrix<int> A(N, N);  //matrix of equation
    
    for( int i = 0; i < N; i++ )
    {
        x(i) = rand() % 3;
    }

    for( int i = 0; i < N; i++ )
    {
        for( int j = 0; j < N; j++ )
        {
            A(i, j) = rand() % 3;
        }
    }
    
    A = prod(trans(A), A);
    f = prod(A, x);
    
    write_data_to_file("data.txt", N, E, A, f);
    write_answer_to_file("answers.txt", x);
    
    return 0;
}
