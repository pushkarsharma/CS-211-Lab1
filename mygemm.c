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
    int i, j, k, l;
    for (i = 0; i < n; i += 3)
    {
        for (j = 0; j < n; j += 3)
        {
            register double c_00 = C[i * n + j];
            register double c_01 = C[i * n + (j + 1)];
            register double c_02 = C[i * n + (j + 2)];

            register double c_10 = C[(i + 1) * n + j];
            register double c_11 = C[(i + 1) * n + (j + 1)];
            register double c_12 = C[(i + 1) * n + (j + 2)];

            register double c_20 = C[(i + 2) * n + j];
            register double c_21 = C[(i + 2) * n + (j + 1)];
            register double c_22 = C[(i + 2) * n + (j + 2)];

            for (k = 0; k < n; k += 3)
            {
                for (l = 0; l < 3; l++)
                {
                    register double a_0l = A[i * n + k + l];
                    register double a_1l = A[(i + 1) * n + k + l];
                    register double a_2l = A[(i + 2) * n + k + l];

                    register double b_l0 = B[(k + l) * n + j];
                    register double b_l1 = B[(k + l) * n + j + 1];
                    register double b_l2 = B[(k + l) * n + j + 2];

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

            C[i * n + j] = c_00;
            C[i * n + (j + 1)] = c_01;
            C[i * n + (j + 2)] = c_02;

            C[(i + 1) * n + j] = c_10;
            C[(i + 1) * n + (j + 1)] = c_11;
            C[(i + 1) * n + (j + 2)] = c_12;

            C[(i + 2) * n + j] = c_20;
            C[(i + 2) * n + (j + 1)] = c_21;
            C[(i + 2) * n + (j + 2)] = c_22;
        }
    }
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
    int i, j, k;
    for (k = 0; k < n; k += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (j = 0; j < n; j += b)
            {
                int i1, j1, k1;
                for (k1 = k; k1 < k + b; k1++)
                {
                    for (i1 = i; i1 < i + b; i1++)
                    {
                        register double r = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b; j1++)
                        {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }
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
    int i, j, k;
    for (i = 0; i < n; i += b)
    {
        for (k = 0; k < n; k += b)
        {
            for (j = 0; j < n; j += b)
            {
                int i1, j1, k1;
                for (i1 = i; i1 < i + b; i1++)
                {
                    for (k1 = k; k1 < k + b; k1++)
                    {
                        register double r = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b; j1++)
                        {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }
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
    int i, j, k;
    for (j = 0; j < n; j += b)
    {
        for (k = 0; k < n; k += b)
        {
            for (i = 0; i < n; i += b)
            {
                int i1, j1, k1;
                for (j1 = j; j1 < j + b; j1++)
                {
                    for (k1 = k; k1 < k + b; k1++)
                    {
                        register double r = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                        }
                    }
                }
            }
        }
    }
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
    int i, j, k;
    for (k = 0; k < n; k += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (i = 0; i < n; i += b)
            {
                int i1, j1, k1;
                for (k1 = k; k1 < k + b; k1++)
                {
                    for (j1 = j; j1 < j + b; j1++)
                    {
                        register double r = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                        }
                    }
                }
            }
        }
    }
}
//Cache Reuse part 3 End

//Cache Reuse part 4
void optimal(const double *A, const double *B, double *C, const int n, const int b)
{
}
