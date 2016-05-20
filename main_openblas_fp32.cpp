#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctime>
#include <iostream>
#include <cstdlib>
extern "C"
{
	#include <cblas.h> // OpenBLAS
}

#define re (0)
#define im (1)
#define rounds (1000000)
#define A_R (3) // The number of rows and columns in each matrix.
#define A_C (3) // This makes the arguments in OpenBLAS easier to understand.
#define B_R (3) // Must be same as A_C
#define B_C (3)
#define V_R (3)

/*
 *   multiply complex 3x3 matrix
 *        C = A X B
 */

void multiply_complex_matrix(float A[3][3][2], float B[3][3][2], float C[3][3][2])
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
void multiply_complex_matvec(float A[3][3][2],float V[3][2],float W[3][2])
{
  int i;
  for(i=0;i<3;i++) {
    W[i][re] = A[i][0][re]*V[0][re]-A[i][0][im]*V[0][im]+
               A[i][1][re]*V[1][re]-A[i][1][im]*V[1][im]+
               A[i][2][re]*V[2][re]-A[i][2][im]*V[2][im];

    W[i][im] = A[i][0][re]*V[0][im]+A[i][0][im]*V[0][re]+
               A[i][1][re]*V[1][im]+A[i][1][im]*V[1][re]+
               A[i][2][re]*V[2][im]+A[i][2][im]*V[2][re];
  }
}

float fRand(float fMin, float fMax)
{
    float f = (float)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void populateRandMatFlat(float A[A_R*A_C*2]) {
    for(int i=0; i<A_R; i++)
        for(int j=0; j<A_C; j++)
            for(int k=0; k<2; k++)
                A[i*A_C*2+j*2+k] = fRand(0,1000);
}

void populateRandVecFlat(float V[V_R*2]) {
    for(int i=0; i<3; i++)
        for(int j=0; j<2; j++)
            V[i*2+j] = fRand(0,1000);
}

void populateRandMat(float A[A_R][A_C][2]) {
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            for(int k=0; k<2; k++)
                A[i][j][k] = fRand(0,1000);
}

void populateRandVec(float V[V_R][2]) {
    for(int i=0; i<3; i++)
        for(int j=0; j<2; j++)
            V[i][j] = fRand(0,1000);
}

int main()
{
    std::srand(0);
    float *A, *B, *C;
    float alpha[2] = {1.0, 0.0};
	float beta[2] = {1.0, 0.0};

    A = (float *) malloc( A_R*A_C*2*sizeof( float ) );
    B = (float *) malloc( B_R*B_C*2*sizeof( float ) );
    C = (float *) malloc( A_R*B_C*2*sizeof( float ) );
    if (A == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      free(A);
      free(B);
      free(C);
      return 1;
    }

	long elapsed_mm = 0;
	clock_t begin, end;
   	begin = clock();
    for(int r=0; r<rounds; r++)
    {
		populateRandMatFlat(A);
		populateRandMatFlat(B);

    	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				A_R, B_C, A_C, alpha, A, B_R, B, B_C, beta, C, B_C);
	}
   	end = clock();
   	elapsed_mm += end-begin;

    //printf ("\n Computations completed.\n");
    printf ("CPU ops: %d\n", elapsed_mm);
    printf ("Time   : %e sec\n", float(elapsed_mm)/CLOCKS_PER_SEC);

    //printf ("\n Deallocating memory \n\n");
    free(A);
    free(B);
    free(C);

    return 0;
}
