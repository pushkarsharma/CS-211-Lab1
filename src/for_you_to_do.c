#include "../include/for_you_to_do.h"

int get_block_size()
{
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
    double *tempv = (double *)malloc(sizeof(double) * n);
    for (i = 0; i < n; i++)
    {
        maxIndex = i;
        max = fabs(A[i * n + i]);

        int t;
        for (t = i + 1; t < n; t++)
        {
            if (fabs(A[t * n + i]) > max)
            {
                maxIndex = t;
                max = fabs(A[t * n + i]);
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
                int temp = ipiv[i];
                ipiv[i] = ipiv[maxIndex];
                ipiv[maxIndex] = temp;
                memcpy(tempv, A + i * n, n * sizeof(double));
                memcpy(A + i * n, A + maxIndex * n, n * sizeof(double));
                memcpy(A + maxIndex * n, tempv, n * sizeof(double));
            }
        }

        for (t = i + 1; t < n; t++)
        {
            A[t * n + i] = A[t * n + i] / A[i * n + i];
            int k;
            for (k = i + 1; k < n; k++)
            {
                A[t * n + k] -= A[t * n + i] * A[i * n + k];
            }
        }
    }
    free(tempv);
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
    double *y = (double *)malloc(n * sizeof(double));
    int i, j;
    double sum;
    if (UPLO == 'L')
    {
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++)
        {
            sum = 0.0;
            for (j = 0; j < i; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    else if (UPLO == 'U')
    {
        y[n - 1] = B[n - 1] / A[(n - 1) * n + n - 1];
        for (i = n - 2; i >= 0; i--)
        {
            sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = (B[i] - sum) / A[i * n + i];
        }
    }

    memcpy(B, y, sizeof(double) * n);
    free(y);
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
    int l;
    for (new_i= i; new_i< n; new_i+= 3)
    {
        for (new_j = j; new_j < n; new_j += 3)
        {
            register double c_00 = C[new_i* n + new_j];
            register double c_01 = C[new_i* n + (new_j + 1)];
            register double c_02 = C[new_i* n + (new_j + 2)];

            register double c_10 = C[(new_i+ 1) * n + new_j];
            register double c_11 = C[(new_i+ 1) * n + (new_j + 1)];
            register double c_12 = C[(new_i+ 1) * n + (new_j + 2)];

            register double c_20 = C[(new_i+ 2) * n + new_j];
            register double c_21 = C[(new_i+ 2) * n + (new_j + 1)];
            register double c_22 = C[(new_i+ 2) * n + (new_j + 2)];

            for (new_k = k; new_k < n; new_k += 3)
            {
                for (l = 0; l < 3; l++)
                {
                    register double a_0l = A[new_i* n + new_k + l];
                    register double a_1l = A[(new_i+ 1) * n + new_k + l];
                    register double a_2l = A[(new_i+ 2) * n + new_k + l];

                    register double b_l0 = B[(new_k + l) * n + new_j];
                    register double b_l1 = B[(new_k + l) * n + new_j + 1];
                    register double b_l2 = B[(new_k + l) * n + new_j + 2];

                    c_00 += a_0l * b_l0;
                    c_01 += a_0l * b_l1;
                    c_02 += a_0l * b_l2;

                    c_10 += a_1l * b_l0;
                    c_11 += a_1l * b_l1;
                    c_12 += a_1l * b_l2;

                    c_20 += a_2l * b_l0;
                    c_21 += a_2l * b_l1;
                    c_22 += a_2l * b_l2;
                }
            }

            C[new_i* n + new_j] = c_00;
            C[new_i* n + (new_j + 1)] = c_01;
            C[new_i* n + (new_j + 2)] = c_02;

            C[(new_i+ 1) * n + new_j] = c_10;
            C[(new_i+ 1) * n + (new_j + 1)] = c_11;
            C[(new_i+ 1) * n + (new_j + 2)] = c_12;

            C[(new_i+ 2) * n + new_j] = c_20;
            C[(new_i+ 2) * n + (new_j + 1)] = c_21;
            C[(new_i+ 2) * n + (new_j + 2)] = c_22;
        }
    }
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
    int ib, i, j, k, maxIndex;
    double max, sum;
    double *tempv = (double*) malloc(sizeof(double) * n);

    for (ib = 0; ib < n; ib += b)
    {
        for (i = ib; i < ib+b && i < n; i++)
        {
            maxIndex = i;
            max = fabs(A[i*n + i]);
            
            int j;
            for (j = i+1; j < n; j++)
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
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[maxIndex];
                    ipiv[maxIndex] = temp;

                    memcpy(tempv, A + i*n, n * sizeof(double));
                    memcpy(A + i*n, A + maxIndex*n, n * sizeof(double));
                    memcpy(A + maxIndex*n, tempv, n * sizeof(double));
                }
            }

            for (j = i+1; j < n; j++)
            {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                int k;
                for (k = i+1; k < ib+b && k < n; k++)
                {
                    A[j*n + k] -= A[j*n +i] * A[i*n + k];
                }
            }
        }

        for (i = ib; i < ib+b && i < n; i++)
        {
            for (j = ib+b; j < n; j++)
            {
                sum = 0;
                for (k = ib; k < i; k++)
                {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        for (i = ib+b; i < n; i += b)
        {
            for (j = ib+b; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}
