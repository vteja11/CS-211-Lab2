#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, maxIndex;
    double max;
    double *temprow = (double*) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++)
    {
        // pivoting
        maxIndex = i;
        max = fabs(A[i*n + i]);
        
        int j=i+1;
        for (j; j < n; j++)
        {
            if (fabs(A[j*n + i]) > max)
            {
                maxIndex = j;
                max = fabs(A[j*n + i]);
            }
        }
        if (max == 0)
        {
            printf("LU factorization failed: coefficient matrix is singular.\n");
            return -1;
        }
        else
        {
            if (maxIndex != i)
            {
                // save pivoting information
                int temp = ipiv[i];
                ipiv[i] = ipiv[maxIndex];
                ipiv[maxIndex] = temp;
                // swap rows
                memcpy(temprow, A + i*n, n * sizeof(double));
                memcpy(A + i*n, A + maxIndex*n, n * sizeof(double));
                memcpy(A + maxIndex*n, temprow, n * sizeof(double));
            }
        }

        // factorization
        j=i+1;
        for (j; j < n; j++)
        {
            A[j*n + i] = A[j*n + i] / A[i*n + i];
            int k;
            for (k = i+1; k < n; k++)
            {
                A[j*n + k] -= A[j*n +i] * A[i*n + k];
            }
        }
    }
    free(temprow);
    return 0;
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double *y = (double*) malloc(n * sizeof(double));
    int i, j;
    double sum;
    if (UPLO == 'L')
    {
        y[0] = B[PVT[0]];
        i=1;
        for (i; i < n; i++)
        {
            sum = 0.0;
            j=0;
            for (j; j < i; j++)
            {
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[PVT[i]] - sum;
        }
    }
    else if (UPLO == 'U')
    {
        y[n - 1] = B[n - 1] / A[(n-1)*n + n-1];
        i=n-2;
        for (i; i >= 0; i--)
        {
            sum = 0;
            j=i+1;
            for (j ; j < n; j++)
            {
                sum += y[j] * A[i*n + j];
            }
            y[i] = (B[i] - sum) / A[i*n + i];
        }
    }

    memcpy(B, y, sizeof(double) * n);
    free(y);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    return 0;
}

