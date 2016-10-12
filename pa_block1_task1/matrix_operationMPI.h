#ifndef MATRIX_OPERATION_MPI
#define MATRIX_OPERATION_MPI

#include <stdbool.h>
#include "matrix_allocation.h"
double** decomposite_ME(double** matrix, size_t columns, size_t rows, size_t n, size_t m, int size, int rank,
                        size_t* new_rows, size_t* new_columns);
double** decomposite_THEM(double** matrix, size_t columns, size_t rows, size_t n, size_t m, int size, int rank,
                        size_t* new_rows, size_t* new_columns);
double** transfer_data(double** matrix, size_t rows, size_t columns, int rank, int size);
void fullfill_data(double** matrix, size_t rows, size_t columns, double** decomposited, size_t dr, size_t dc,
                   int rank, int size);
void Transpose(double*** matrix, size_t *columns, size_t *rows, size_t n, size_t m, int size, int rank);
double** DistributeColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int starter, size_t* result_rows);
void CollectColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int collector,
                    double** result, size_t global_rows, size_t global_cols);

/// implementation
/// in each proc: matrix column and rows
/// n and m are for the whole big (global) matrix
void Transpose(double*** matrix, size_t *columns, size_t *rows, size_t n, size_t m, int size, int rank)
{
    size_t mine_rows, mine_columns;
    double** mine = decomposite_ME(*matrix, *columns, *rows, n, m,size, rank, &mine_rows, &mine_columns);

    size_t them_rows, them_columns;
    double** them = decomposite_THEM(*matrix, *columns, *rows, n, m,size, rank, &them_rows, &them_columns);
    //printf("-%d Entering transfer\n", rank);
    double** result = transfer_data(them, them_rows, them_columns, rank, size);
    //printf("-%d leaving transfer\n", rank);
    fullfill_data(mine, mine_rows, mine_columns, result, them_rows, them_columns, rank, size);

    free_matrix(matrix, *rows, *columns);

    free_matrix(&them, them_rows, them_columns);

    free_matrix(&result, them_rows, them_columns);
    *matrix = mine;
    *rows = mine_rows;
    *columns = mine_columns;
}

/// each proc: returns the transpoed matrix with few fullfiled cells from the data inside of the processor
/// n and m are for the whole big matrix
double** decomposite_ME(double** matrix, size_t columns, size_t rows, size_t n, size_t m, int size, int rank,
                        size_t* new_rows, size_t* new_columns)
{
    *new_rows = m / size + ((rank < (m % size)) ? 1 : 0);
    *new_columns = n;
    double** result = allocate_matrix(*new_rows, n);

    size_t i = 0, j = 0;
    for(; i < rows; ++i)
        for(j = 0; j*size + rank < columns; ++j)
            result[j][i*size+rank] = matrix[i][j*size + rank];

    return result;
}
/// returns matrix, every row represents the proc to whome to send
double** decomposite_THEM(double** matrix, size_t columns, size_t rows, size_t n, size_t m, int size, int rank,
                        size_t* new_rows, size_t* new_columns)
{
    *new_rows = size;
    size_t add_r = (n > m) ? 1 : 0;
    size_t add_c = (m >= n) ? 1 : 0;
    *new_columns = (n/size + add_c) * ((int)(columns/size) + add_r + 1);  // not working
    //if(rank == 0)
    //    *new_columns = 2;
    //printf("\n----------%d--------%d,%d, %d, %d, (%d,%d), \n",rank, *new_columns , rows, columns, n, add_r, add_c);
    double** result = allocate_matrix(size, *new_columns);
    size_t proc = 0,i = 0, j, count=0;

    for(; proc < size; ++proc)
    {
        if(proc == rank)
            continue;
        count = 0;
        for(j = proc; j < columns ; j+=size)
            for(i = 0; i < rows; ++i) // interchange loops nahh
                result[proc][count++] = matrix[i][j];
    }
    return result;
}

double** transfer_data(double** matrix, size_t rows, size_t columns, int rank, int size)
{
    //printf("=--%d===\n", rank);
    double** temporary = allocate_matrix(rows, columns);
    int ierr, seq, sign, bound = (size-1)/2, inc=((size-1)%2);
    MPI_Status status;
    size_t k;
    seq = 0;
    sign = 1;
    for(k = 1; k < size; ++k)
    {
        seq += sign * 2;
        if(seq > bound*2)
        {
            if(!inc)
                seq -= 2;
            sign = -1;
        }
        //printf("---%d %d %d, k=%d\n", rank, seq, bound, k);

        //printf("%d, %d, %d, %d, %d \n", rank, seq, k, sign, bound);

        if(rank % seq < seq/2)
        {
            //printf("-k=%d rank=%d seq=%d \n", k, rank, seq);
            ierr = MPI_Send(matrix[(rank + k) % size], columns, MPI_DOUBLE, (rank + k) % size, 0, MPI_COMM_WORLD);
            ierr = MPI_Recv(temporary[(rank - k+size) % size], columns, MPI_DOUBLE, (rank - k+size) % size, 0, MPI_COMM_WORLD, &status);
        }
        else
        {
            ierr = MPI_Recv(temporary[(rank - k+size) % size], columns, MPI_DOUBLE, (rank - k+size) % size, 0, MPI_COMM_WORLD, &status);
            ierr = MPI_Send(matrix[(rank + k) % size], columns, MPI_DOUBLE, (rank + k) % size, 0, MPI_COMM_WORLD);
        }
    }
    return temporary;
}

void fullfill_data(double** matrix, size_t rows, size_t columns, double** decomposited, size_t dr, size_t dc,
                   int rank, int size)
{
    size_t k = 0,i,j,count;
    for(;k < size; ++k)
    {
        if(k == rank)
            continue;
        for(i = 0, count = 0; i < rows; ++i)
        {
            for(j = k; j < columns; j += size)
            {
                matrix[i][j] = decomposited[k][count++];
            }
        }
    }
}


#if 1
/// it distributes a matrix  1 2 3      p0 1 2 3  p1  4 5 6
/// for ex. 2 procs          4 5 6  ->     7 8 9
/// by rows for              7 8 9
/// for faster allocation
/// packed data
double** DistributeColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int starter, size_t* result_rows)
{
    MPI_Status status;
    double** result = NULL;
    int ierr, i, proc;
    *result_rows = (int)(rows/size) + ((rank < (rows % size)) ? 1 : 0);
    double* communication = NULL;

    result = allocate_matrix(*result_rows, cols);
    if(result == NULL)
        return NULL;

    if(rank == starter)
    {
        communication = malloc(sizeof(double) * ((int)(rows/size) + 1 ) * cols);
        for(proc = 0; proc < size; ++proc)
        {
            if(proc == rank)
            {
                for(i = rank; i < rows; i += size)
                    memcpy(result[i/size], matrix[i], sizeof(double) * cols);
                continue;
            }
            for(i = proc; i < rows; i+=size)
                memcpy(&communication[(i/size)*cols], matrix[i], sizeof(double) * cols);
            ierr = MPI_Send(communication, ((i-1)/size)*cols, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        communication = malloc(sizeof(double) * *result_rows * cols);
        ierr = MPI_Recv(communication, *result_rows * cols, MPI_DOUBLE, starter, 0, MPI_COMM_WORLD, &status);
        for(i = 0; i < *result_rows; ++i)
            memcpy(result[i], communication + i*cols, sizeof(double) * cols);
    }
    free(communication);
    return result;
}

/// here rows,cols means the rows,cols in each processor not the global rows,cols of the matrix
void CollectColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int collector,
                    double** result, size_t global_rows, size_t global_cols)
{
    MPI_Status status;
    int ierr, i, proc, comm_size;
    double* communication = NULL;

    if(rank == collector)
    {
        communication = malloc(sizeof(double) * (global_rows/size + 1) * cols);
        for(proc = 0; proc < size; ++proc)
        {
            if(proc == rank)
            {
                for(i = rank; i < global_rows; i += size)
                    memcpy(result[i], matrix[i/size], sizeof(double) * cols);

                continue;
            }
            comm_size = global_rows/size + ((proc < global_rows%size) ? 1 : 0);

            ierr = MPI_Recv(communication, comm_size*cols, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);

            for(i = 0; i < comm_size; ++i)
                memcpy(result[proc+i*size], communication + i*cols, sizeof(double) * cols);
        }
        free(communication);
    }
    else
    {
        communication = malloc(sizeof(double) * rows * cols);
        for(i = 0; i < rows; ++i)
            memcpy(communication + i*cols, matrix[i], sizeof(double) * cols);

        ierr = MPI_Send(communication, rows * cols, MPI_DOUBLE, collector, 0, MPI_COMM_WORLD);
        free(communication);
    }
}
#else
/// it distributes a matrix  1 2 3      p0 1 2 3  p1  4 5 6
/// for ex. 2 procs          4 5 6  ->     7 8 9
/// by rows for              7 8 9
/// for faster allocation
double** DistributeColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int starter, size_t* result_rows)
{
    MPI_Status status;
    double** result = NULL;
    int ierr, i;
    *result_rows = (int)(rows/size) + ((rank < (rows % size)) ? 1 : 0);

    result = allocate_matrix(*result_rows, cols);
    if(result == NULL)
        return NULL;

    if(rank == starter)
    {
        for(i = 0;i < rows; ++i)
        {
            if(i%size == starter)
            {
                memcpy(result[(int)i/size], matrix[i], sizeof(double) * cols);
                continue;
            }
            ierr = MPI_Send(matrix[i], cols, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        for(i = 0; i < *result_rows; ++i)
            ierr = MPI_Recv(result[i], cols, MPI_DOUBLE, starter, 0, MPI_COMM_WORLD, &status);
    }
    return result;
}

/// here rows,cols means the rows,cols in each processor not the global rows,cols of the matrix
void CollectColumns(double** matrix, size_t rows, size_t cols, int size, int rank, int collector,
                    double** result, size_t global_rows, size_t global_cols)
{
    MPI_Status status;
    int ierr, i;

    if(rank == collector)
    {
        for(i = 0;i < global_rows; ++i)
        {
            if(i%size == collector)
            {
                memcpy(result[i], matrix[(int)(i/size)], sizeof(double) * cols);
                continue;
            }
            ierr = MPI_Recv(result[i], cols, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        for(i = 0; i < rows; ++i)
            ierr = MPI_Send(matrix[i], cols, MPI_DOUBLE, collector, 0, MPI_COMM_WORLD);
    }
}
#endif
#endif
