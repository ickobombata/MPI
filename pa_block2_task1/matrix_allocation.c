#include "matrix_allocation.h"

void copy_Matrix1D(struct Matrix1D* mat, struct Matrix1D* source)
{
    int i, mn = mat->m * mat->n;
    for(i = 0 ; i < mn; ++i)
        mat->data[i] = source->data[i];
}

bool allocate_Matrix1D(size_t rows, size_t columns, struct Matrix1D* matrix)
{
    matrix->n = rows;
    matrix->m = columns;
    matrix->data = (double*)malloc(sizeof(double) * rows * columns);
    if(matrix->data == NULL)
        return false;

    return true;
}

void free_Matrix1D(struct Matrix1D* matrix)
{
    free(matrix->data);
}

double* get_Matrix1D(struct Matrix1D* matrix, size_t i, size_t j)
{
    return &matrix->data[i * matrix->m + j];
}


void print_Matrix1D(struct Matrix1D* matrix)
{
    int i,j;
    for(i = 0; i < matrix->n; ++i)
    {
        for(j = 0; j < matrix->m; ++j)
            printf("%6.2lf, ", matrix->data[i*matrix->m + j]);
        printf("\n");
    }

//    for(i = 0; i < 12; ++i)
//    {
//        for(j = 0; j < 12; ++j)
//            printf("%6.2lf, ", matrix->data[i*matrix->m + j]);
//        printf("\n");
//    }
}

double** allocate_matrix(size_t rows, size_t columns)
{
    double** matrix = NULL;
    if((matrix = malloc(sizeof(double*) * rows)) == NULL)
        return NULL;
    int i = 0;
    for(;i < rows; ++i)
    {
        if((matrix[i] = malloc(sizeof(double) * columns)) == NULL)
        {
            while(--i >= 0)
            {
                free(matrix[i]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

void free_matrix(double*** matrix, size_t rows, size_t columns)
{
    int i = 0;
    for(; i < rows; ++i)
        free((*matrix)[i]);
    free(*matrix);
    *matrix = NULL;
}

void print_matrix(double** matrix, size_t rows, size_t columns)
{
    int i =0, j=0;
    for(;i < rows; ++i)
    {
        for(;j< columns; ++j)
            printf("%2.4f,", matrix[i][j]);
        printf("\n");
        j = 0;
    }
}

void print_vector(double* vec, size_t elements)
{
    int i = 0;
    for(; i < elements; ++i)
        printf("%2.4lf, ", vec[i]);
    printf("\n");
}

double* matrix_to_vector(double** matrix, size_t rows, size_t columns)
{
    double* result = (double*)malloc(sizeof(double) * rows * columns);
    size_t i = 0;
    for(i = 0 ; i < rows; ++i)
        memcpy(result + i * columns, matrix[i], columns * sizeof(double));
    return result;
}

double** vector_to_matrix(double* arr, size_t rows, size_t columns)
{
    double** result = allocate_matrix(rows, columns);
    size_t i = 0;
    for(i = 0 ; i < rows; ++i)
        memcpy(result[i], arr + i * columns, columns * sizeof(double));
    return result;
}

