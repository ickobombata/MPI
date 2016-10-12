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
            printf("%2.5f ,", matrix[i][j]);
        printf("\n");
        j = 0;
    }
}

#endif
