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
    int i1, j1, k1, l;
    for (i1 = i; i1 < (i + b > n ? n : (i + b)); i1 += 3)
    {
        for (j1 = j; j1 < (j + b > n ? n : (j + b)); j1 += 3)
        {
            register double c00 = C[i1 * n + j1];
            register double c01 = C[i1 * n + (j1 + 1)];
            register double c02 = C[i1 * n + (j1 + 2)];

            register double c10 = C[(i1 + 1) * n + j1];
            register double c11 = C[(i1 + 1) * n + (j1 + 1)];
            register double c12 = C[(i1 + 1) * n + (j1 + 2)];

            register double c20 = C[(i1 + 2) * n + j1];
            register double c21 = C[(i1 + 2) * n + (j1 + 1)];
            register double c22 = C[(i1 + 2) * n + (j1 + 2)];

            for (k1 = k; k1 < ((k + b) > n ? n : (k + b)); k1 += 3)
            {
                for (l = 0; l < 3; l++)
                {
                    //accessing A column wise
                    register double a_0l = A[i1 * n + k1 + l];       //row1
                    register double a_1l = A[(i1 + 1) * n + k1 + l]; //row2
                    register double a_2l = A[(i1 + 2) * n + k1 + l]; //row3

                    //accessing B row wise
                    register double b_l0 = B[(k1 + l) * n + j1];     //col1
                    register double b_l1 = B[(k1 + l) * n + j1 + 1]; //col2
                    register double b_l2 = B[(k1 + l) * n + j1 + 2]; //col3

                    c00 -= a_0l * b_l0;
                    c01 -= a_0l * b_l1;
                    c02 -= a_0l * b_l2;

                    c10 -= a_1l * b_l0;
                    c11 -= a_1l * b_l1;
                    c12 -= a_1l * b_l2;

                    c20 -= a_2l * b_l0;
                    c21 -= a_2l * b_l1;
                    c22 -= a_2l * b_l2;
                }
            }

            C[i1 * n + j1] = c00;
            C[i1 * n + (j1 + 1)] = c01;
            C[i1 * n + (j1 + 2)] = c02;

            C[(i1 + 1) * n + j1] = c10;
            C[(i1 + 1) * n + (j1 + 1)] = c11;
            C[(i1 + 1) * n + (j1 + 2)] = c12;

            C[(i1 + 2) * n + j1] = c20;
            C[(i1 + 2) * n + (j1 + 1)] = c21;
            C[(i1 + 2) * n + (j1 + 2)] = c22;
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
    int i, j, k, t, maxind, ib;
    double max;
    double *tempv = (double *)malloc(sizeof(double) * n);
    for (ib = 0; ib < n; ib += b)
    {
        for (i = ib; i < ib + b && i < n; i++)
        {
            maxind = i;
            max = fabs(A[i * n + i]);

            for (t = i + 1; t < n; t++)
            {
                if (fabs(A[t * n + i]) > max)
                {
                    maxind = t;
                    max = fabs(A[t * n + i]);
                }
            }

            if (max == 0)
            {
                printf("LUfactoration failed: coefficient matrix is singular\n");
                return -1;
            }
            else
            {
                if (maxind != i)
                {
                    int temps = ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = temps;

                    memcpy(tempv, A + i * n, n * sizeof(double));
                    memcpy(A + i * n, A + maxind * n, n * sizeof(double));
                    memcpy(A + maxind * n, tempv, n * sizeof(double));
                }

                for (j = i + 1; j < n; j++)
                {
                    A[j * n + i] = A[j * n + i] / A[i * n + i];
                    for (k = i + 1; k < ib + b && k < n; k++)
                        A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
                }
            }
        }
        for (i = ib; i < ib + b && i < n; i++)
        {
            for (j = ib + b; j < n; j++)
            {
                double sum = 0.0;
                for (k = ib; k < i; k++)
                {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] -= sum;
            }
        }
        for (i = ib + b; i < n; i += b)
        {
            for (j = ib + b; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}
