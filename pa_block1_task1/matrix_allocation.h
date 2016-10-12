#ifndef MATRIX_ALLOCATION
#define MATRIX_ALLOCATION

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

double** vector_to_matrix(double* vector, size_t rows, size_t columns)
{
    double** result = allocate_matrix(rows, columns);
    size_t i = 0;
    for(i = 0 ; i < rows; ++i)
        memcpy(result[i], vector + i * columns, columns * sizeof(double));
    return result;
}


#endif
