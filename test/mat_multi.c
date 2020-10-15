#include <stdio.h>
#include <stdbool.h>

void print(int *matrix, int n)
{
    int pos;
    for (pos = 0; pos < n * n; pos++)
    {
        printf("%d\t", *(matrix + pos));
        if ((pos + 1) % n == 0)
        {
            printf("\n");
        }
    }
    printf("----------------------------\n");
}

void ijk(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register int sum = 0;
            for (k = 0; k < n; k++)
            {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
    print(C, 4);
}

void jik(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            register int sum = 0;
            for (k = 0; k < n; k++)
            {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
    print(C, 4);
}

void kij(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            register int r = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
    print(C, 4);
}

void ikj(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            register int r = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
    print(C, 4);
}

void jki(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            register int r = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }
    print(C, 4);
}

void kji(const int *A, const int *B, int *C, const int n)
{
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            register int r = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * r;
            }
        }
    }
    print(C, 4);
}

void bijk(const int *A, const int *B, int *C, const int n, const int b)
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

    print(C, 4);
}

void bjik(const int *A, const int *B, int *C, const int n, const int b)
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
    print(C, 4);
}

int main()
{
    int A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    print(A, 4);
    int B[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    int C[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int D[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int E[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int F[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int G[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int H[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // int A[] = {1, 2, 3, 4};
    // int B[] = {1, 2, 3, 4};
    // int C[] = {0, 0, 0, 0};

    // print(A, 2);

    ijk(A, B, C, 4);
    // jik(A, B, D, 4);
    // kij(A, B, E, 4);
    // ikj(A, B, F, 4);
    // jki(A, B, G, 4);
    // kji(A, B, H, 4);

    bjik(A, B, D, 4, 2);

    // bool all = true;
    // int pos;
    // for (pos = 0; pos < 16; pos++)
    // {
    //     if (*(C + pos) != *(D + pos))
    //     {
    //         all = false;
    //         break;
    //     }
    // }
    // if (all)
    // {
    //     printf("true");
    // }
    // else
    // {
    //     printf("false");
    // }
    // printf("\n");

    return 0;
}