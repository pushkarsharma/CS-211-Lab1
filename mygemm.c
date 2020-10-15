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
    int i, j, k;

    for (i = 0; i < n; i += 2)
    {
        for (j = 0; j < n; j += 2)
        {
            register double C_00 = C[i * n + j];
            register double C_01 = C[i * n + (j + 1)];
            register double C_10 = C[(i + 1) * n + j];
            register double C_11 = C[(i + 1) * n + (j + 1)];

            for (k = 0; k < n; k += 2)
            {
                register double A_00 = A[i * n + k];
                register double A_01 = A[i * n + (k + 1)];
                register double A_10 = A[(i + 1) * n + k];
                register double A_11 = A[(i + 1) * n + (k + 1)];

                register double B_00 = B[k * n + j];
                register double B_01 = B[k * n + (j + 1)];
                register double B_10 = B[(k + 1) * n + j];
                register double B_11 = B[(k + 1) * n + (j + 1)];

                C_00 += A_00 * B_00 + A_01 * B_10;
                C_01 += A_00 * B_01 + A_01 * B_11;
                C_10 += A_10 * B_00 + A_11 * B_10;
                C_11 += A_10 * B_01 + A_11 * B_11;
            }
            C[i * n + j] = C_00;
            C[i * n + (j + 1)] = C_01;
            C[(i + 1) * n + j] = C_10;
            C[(i + 1) * n + (j + 1)] = C_11;
        }
    }
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
    int i, j, k;
    for (j = 0; j < n; j += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (k = 0; k < n; k += b)
            {
                int i1, j1, k1;
                for (j1 = j; j1 < j + b; j1++)
                {
                    for (i1 = i; i1 < i + b; i1++)
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
