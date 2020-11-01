#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */

    // for testing purpose
    //int a[] ={99,129,132,150,177,270};
    //return a[i];
    return 132;
  
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
    
    int i, max_i;
    double max;
    double *row_tmp = (double*) malloc(sizeof(double) * n);
    
    //
    i=0;
    for (i; i < n; i++)
    {
        // we are pivoting in this step
        max_i = i;
        max = fabs(A[i*n + i]);
        
        int j=i+1;
        for (j; j < n; j++)
        {
            if (fabs(A[j*n + i]) > max)
            {
                max_i = j;
                max = fabs(A[j*n + i]);
            }
        }

        //check for singular matrix case after pivoting
        if (max == 0)
        {     
            return -1;
        }
        else
        {
            if (max_i != i)
            {
                // pivoting information is strored in this step
                int temp = ipiv[i];
                ipiv[i] = ipiv[max_i];
                ipiv[max_i] = temp;

                // swapiing the rows
                memcpy(row_tmp, A + i*n, n * sizeof(double));
                memcpy(A + i*n, A + max_i*n, n * sizeof(double));
                memcpy(A + max_i*n, row_tmp, n * sizeof(double));
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
    free(row_tmp);
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

    //if lower triangular matrix
    if (UPLO == 'L')
    {
        y[0] = B[ipiv[0]];
        i=1;
        for (i; i < n; i++)
        {
            sum = 0.0;
            j=0;
            for (j; j < i; j++)
            {
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    //if upper triangular matrix
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

    //code from lab1
    int i_tmp = i, j_tmp = j, k_tmp = k;

    //boundary for i
    int n1 = i + b > n ? n : i + b;
    //boundary for j
    int n2 = j + b > n ? n : j + b;
    //boundary for k
    int n3 = k + b > n ? n : k + b;

    for (i_tmp = i; i_tmp < n1; i_tmp += 3)
    {
        for (j_tmp = j; j_tmp < n2; j_tmp += 3)
        {
            int t = i_tmp * n + j_tmp;
            int t2 = t + n;
            int t3 = t2 + n;

            //registers for C
            register double c1 = C[t];
            register double c2 = C[t + 1];
            register double c3 = C[t + 2];
            register double c4 = C[t2];
            register double c5 = C[t2 + 1];
            register double c6 = C[t2 + 2];
            register double c7 = C[t3];
            register double c8 = C[t3 + 1];
            register double c9 = C[t3 + 2];

            for (k_tmp = k; k_tmp < n3; k_tmp += 3)
            {
                int l;
                for (l = 0; l < 3; l++)
                {
                    int ta = i_tmp * n + k_tmp + l;
                    int t2a = ta + n;
                    int t3a = t2a + n;
                    int tb = k_tmp * n + j_tmp + l * n;

                    // registers for A 
                    register double a0 = A[ta];
                    register double a1 = A[t2a];
                    register double a2 = A[t3a];
                    // registers for B
                    register double b0 = B[tb];
                    register double b1 = B[tb + 1];
                    register double b2 = B[tb + 2];

                    c1 -= a0 * b0;
                    c2 -= a0 * b1;
                    c3 -= a0 * b2;
                    c4 -= a1 * b0;
                    c5 -= a1 * b1;
                    c6 -= a1 * b2;
                    c7 -= a2 * b0;
                    c8 -= a2 * b1;
                    c9 -= a2 * b2;
                }
            }
            C[t] = c1;
            C[t + 1] = c01;
            C[t + 2] = c3;
            C[t2] = c4;
            C[t2 + 1] = c5;
            C[t2 + 2] = c6;
            C[t3] = c7;
            C[t3 + 1] = c8;
            C[t3 + 2] = c9;
        }
    }
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
    int i_block,  max_i, i, j, k;
    
    double *row_tmp = (double*) malloc(sizeof(double) * n);
    double max, sum;

    for (i_block = 0; i_block < n; i_block += b)
    {
        for (i = i_block; i < i_block+b && i < n; i++)
        {
            // we are pivoting in this step
            max_i = i;
            max = fabs(A[i*n + i]);
            
            int j;
            for (j = i+1; j < n; j++)
            {
                if (fabs(A[j*n + i]) > max)
                {
                    max_i = j;
                    max = fabs(A[j*n + i]);
                }
            }

            //check for singular matrix 
            if (max == 0)
            {
                return -1;
            }
            else
            {
                if (max_i != i)
                {
                    // pivoting information is stored here
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[max_i];
                    ipiv[max_i] = temp;
                    
                    // swap rows
                    memcpy(row_tmp, A + i*n, n * sizeof(double));
                    memcpy(A + i*n, A + max_i*n, n * sizeof(double));
                    memcpy(A + max_i*n, row_tmp, n * sizeof(double));
                }
            }

            // factorization
            for (j = i+1; j < n; j++)
            {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                int k;
                for (k = i+1; k < i_block+b && k < n; k++)
                {
                    A[j*n + k] -= A[j*n +i] * A[i*n + k];
                }
            }
        }

        // updating  A from A[i,j] as i=i_block:end of block, j=end_of_block+1:n
        for (i = i_block; i < i_block+b && i < n; i++)
        {
            for (j = i_block+b; j < n; j++)
            {
                sum = 0;
                for (k = i_block; k < i; k++)
                {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        // updating A from A[i,j] as i=end of block+1:n, end of block+1:n)
        for (i = i_block+b; i < n; i += b)
        {
            for (j = i_block+b; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, i_block, b);
            }
        }
    }
    return 0;
}

