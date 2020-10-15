#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double r = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n)
{
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n)
{
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double r = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    for (i = 0; i < n; i += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (k = 0; k < n; k += b)
            {
                int i1, j1, k1;
                for (i1 = i; i1 < i + b; i1++)
                {
                    for (j1 = j; j1 < j + b; j1++)
                    {
                        register int r = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b; k1++)
                        {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            register double r = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b)
{
    /*
    int i, j, k;
    for (j = 0; i < n; i++)
    {
        for (i = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                int i1, j1, k1;
                for (j1 = j; j1 < j + b; j1++)
                {
                    for (i1 = j; i1 < i + b; i1++)
                    {
                        register double r = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b; k1++)
                        {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }
*/
}

void kij(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            register double r = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b)
{
}

void ikj(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            register double r = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b)
{
}

void jki(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            register double r = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b)
{
}

void kji(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            register double r = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b)
{
}
//Cache Reuse part 3 End

//Cache Reuse part 4
void optimal(const double *A, const double *B, double *C, const int n, const int b)
{
}
