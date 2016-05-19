#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <complex> // Needed for complex numbers with OpenBLAS (and other BLAS libraries).
extern "C"
{
    #include <cblas.h> // OpenBlas
}

#define re (0)
#define im (1)
#define rounds (5000000) 
#define A_R (3) // The number of rows and columns in each matrix.
#define A_C (3) // This makes the arguments in OpenBLAS easier to understand.
#define B_R (3)
#define B_C (3)

/*
 *   multiply complex 3x3 matrix 
 *        C = A X B
 */
void multiply_complex_matrix(double A[3][3][2],double B[3][3][2],double C[3][3][2])
{
  int i,j,k;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        C[i][j][re] += A[i][k][re]*B[k][j][re]-A[i][k][im]*B[k][j][im];
        C[i][j][im] += A[i][k][im]*B[k][j][re]+A[i][k][re]*B[k][j][im];
      }
    }
  }
}

/*
 *   multiply complex 3x3 matrix and 3 vector
 *        W = A X V
 */
void multiply_complex_matvec(double A[3][3][2],double V[3][2],double W[3][2])
{
  int i;
  for(i=0;i<3;i++) {       
    W[i][re] = A[i][0][re]*V[0][re]-A[i][0][im]*V[0][im]+
               A[i][1][re]*V[1][re]-A[i][1][im]*V[1][im]+
               A[i][2][re]*V[2][re]-A[i][2][im]*V[2][im] ;
    W[i][im] = A[i][0][re]*V[0][im]+A[i][0][im]*V[0][re]+
               A[i][1][re]*V[1][im]+A[i][1][im]*V[1][re]+
               A[i][2][re]*V[2][im]+A[i][2][im]*V[2][re] ;
  }
} 

double fRand(double fMin, double fMax)
{
    double f = (double)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main()
{    
    std::srand(std::time(0));
    double A[3][3][2];
    double B[3][3][2];
    double C[3][3][2];
    double V[3][2];
    double W[3][2];
    double A2[18];
    double B2[18];
    double C2[18];
    double V2[6];
    double W2[6];
    // Generate some random numbers
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<2; k++)
            {
                A[i][j][k] = fRand(0,1000);
                A2[k + j*2 + i*6] = A[i][j][k];
                B[i][j][k] = fRand(0,1000);
                B2[k + j*2 + i*6] = B[i][j][k];
                C[i][j][k] = 0.0;
                C2[k + j*2 + i*6] = 0.0;
            }
        }
    }
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<2; j++)
        {
            V[i][j] = fRand(0,1000);
            V2[i*2+j] = V[i][j];
        }
    }
    
    ////////////// Calculate with mosc3 
    clock_t begin = clock();
    // MxM
    for(int i=0; i<rounds; i++)
    {
        multiply_complex_matrix(A,B,C);
    }
    clock_t end = clock();
    double elapsed_mosc3_mm = double(end-begin)/CLOCKS_PER_SEC;
    
    begin = clock();
    // MxV
    for(int i=0; i<rounds; i++)
    {
        multiply_complex_matvec(A, V, W);
    }
    end = clock();
    double elapsed_mosc3_mv = double(end-begin)/CLOCKS_PER_SEC;
    
    ////////////// Calculate with OpenBlas 
    begin = clock();
    // MxM
    double alpha[2];
    alpha[0] = 1.0;
    alpha[1] = 0.0;
    for(int i=0; i<rounds; i++)
    {
        // Reference: http://www.netlib.org/lapack/explore-html/dc/d17/group__complex16__blas__level3.html#ga4ef748ade85e685b8b2241a7c56dd21c
        // ZGEMM performs 
        // C := alpha*op(A)*op(B)+beta*C
        // where op(X) is either X or transpose of X etc.
        // alpha and beta are scalars.
        // The first three '3' in the argument are the number of: 
        // rows of A2 and C2, 
        // columns of B2 and columns of C2,
        // columns of A2 and rows of B2
        // The '1' declare the first dimension of each matrix as declared in the calling program.
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_R, B_C, B_R, alpha, A2, A_C, B2, B_C, alpha, C2, B_C);
    }
    end = clock();
    double elapsed_openBlas_mm = double(end-begin)/CLOCKS_PER_SEC;
    
    begin = clock();
    // MxV
    double beta[2];
    beta[0] = 0.0;
    beta[1] = 0.0;
    int stepSize = 1;
    for(int i=0; i<rounds; i++)
    {
        // Reference: http://www.netlib.org/lapack/explore-html/dc/dc1/group__complex16__blas__level2.html#gafaeb2abd9fffa7442b938dc384aeaf47
        cblas_zgemv(CblasRowMajor, CblasNoTrans, A_R, A_C, alpha, A2, A_C, V2, stepSize, beta, W2, stepSize);
    }
    end = clock();
    double elapsed_openBlas_mv = double(end-begin)/CLOCKS_PER_SEC;
    
    std::cout << "Evaluated CPU time with Mosc3 (MxM): " << elapsed_mosc3_mm << std::endl;
    std::cout << "Evaluated CPU time with OpenBLAS (MxM): " << elapsed_openBlas_mm << std::endl;
    std::cout << "Evaluated CPU time with Mosc3 (MxV): " << elapsed_mosc3_mv << std::endl;
    std::cout << "Evaluated CPU time with OpenBLAS (MxV): " << elapsed_openBlas_mv << std::endl;
}