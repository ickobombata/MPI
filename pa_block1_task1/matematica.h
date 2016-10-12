#ifndef MATEMATICA
#define MATEMATICA

void swap_double(double* first, double* second)
{
    double temp = *first;
    *first = *second;
    *second = temp;
}

size_t my_rows_from(size_t n, int rank, int size)
{
    return ((n/size) + ((rank < n%size) ? 1 : 0));
}

void swap_int(int* first, int* second)
{
    int temp = *first;
    *first = *second;
    *second = temp;
}

double absolute(double a)
{
    if(a < 0)
        return -a;
    return a;
}

double frobenius_norm(double** matrix, int rows, int cols)
{
    int i, j;
    double result = 0;
    for(i = 0; i < rows; ++i)
        for(j = 0; j < cols; ++j)
            result += matrix[i][j] * matrix[i][j];
    return sqrt(result);
}

double* identity(int n)
{
    double* result = malloc(sizeof(double) * n * n);
    memset(result, 0, sizeof(double) * n * n);
    int i,j= 0;
    for(i = 0; i < n; ++i)
        result[i*n + j++] = 1;
    return result;
}


/// concatenate matrix A and B adds only references to the rows of B in A
/// NOT DEEP COPY
void concatenate_matrix(double*** A, size_t n, size_t m, double** B, size_t q)
{
    double** temp = *A;
    *A = malloc(sizeof(double*) * (n + q));
    int i = 0, j=0;
    for(i = 0; i < n; ++i)
        (*A)[i] = temp[i];
    for(j = 0; j < q; ++j)
        (*A)[i+j] = B[j];

    free(temp);
}

double** construct_permutation_matrix(int* P, int n)
{
    double** result = allocate_matrix(n, n);
    int i, j;
    for(i = 0; i < n; ++i)
        memset(result[i], 0, sizeof(double) * n);

    for(i = 0; i < n; ++i)
        result[i][P[i]] = 1.0;

    return result;
}

#endif

