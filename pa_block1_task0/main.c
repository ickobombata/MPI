#include "/usr/include/mpi/mpi.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <unistd.h>

#include "matrix_allocation.h"
#include "matrix_operationMPI.h"

/*              m   repres      original
matrix = [[1, 4, 7, 10],        1,   2,  3,
       n  [2, 5, 8, 11],        4,   5,  6,
          [3, 6, 9, 12]]   ==   7,   8,  9,
                                10, 11, 12.
*/
struct Matrix
{
    double** matrix;
    size_t rows;
    size_t columns;
};

double** give_temp_matrix(int rows, int cols)
{
    int i, j;
    double** result = allocate_matrix(cols, rows);
    for(i = 0; i < cols; ++i)
        for(j = 0; j < rows; ++j)
            result[i][j] = j + (double)i/1000;

    return result;
}

bool test_distribute(int rows, int columns, int rank, int size)
{
    int i, j;
    double** matrix = NULL;
    if(rank == 0)
        matrix = give_temp_matrix(rows, columns);

    size_t new_rows;
    double** result = DistributeColumns(matrix, columns, rows, size, rank, 0, &new_rows);

    if(rank == 0)
        printf("Distributed %d %d\n", rows, columns);

    free_matrix(&result, new_rows, rows);
    if(rank == 0)
        free_matrix(&matrix, columns ,rows);
    return true;
}

bool test_collect(int rows, int columns, int rank, int size)
{
    int i, j;
    double** matrix = NULL;
    if(rank == 0)
        matrix = give_temp_matrix(rows, columns);

    size_t new_rows;
    double** result = DistributeColumns(matrix, columns, rows, size, rank, 0, &new_rows);

    if(rank == 0)
        free_matrix(&matrix, columns ,rows);

    if(rank == 0)
        matrix = allocate_matrix(columns, rows);

    CollectColumns(result, new_rows, rows, size, rank, 0, matrix, columns, rows);

    if(rank == 0)
        printf("Collected %d %d\n", rows, columns);

    free_matrix(&result, new_rows, rows);
    if(rank == 0)
        free_matrix(&matrix, columns ,rows);

    return true;
}

bool test_decomposite_me(int rows, int cols, int n, int m, int rank, int size)
{
    if(rank != 0)
        return true;

    int i, j;
    double** matrix = allocate_matrix(rows, cols);
    for(i = 0; i < rows; ++i)
        for(j = 0; j < cols; ++j)
            matrix[i][j] = i*cols + j;

    //print_matrix(matrix, rows, cols);

    size_t res_rows, res_cols;
    double** result = decomposite_ME(matrix, cols, rows, n, m, size, rank, &res_rows, &res_cols);
    //printf("\n");
    //print_matrix(result, res_rows, res_cols);

    free_matrix(&matrix, rows, cols);
    free_matrix(&result, res_rows, res_cols);
    return true;
}

bool test_decomposite_them(int rows, int cols, int n, int m, int rank, int size)
{
    if(rank != 0)
        return true;

    int i, j;
    double** matrix = allocate_matrix(rows, cols);
    for(i = 0; i < rows; ++i)
        for(j = 0; j < cols; ++j)
            matrix[i][j] = i*cols + j;

    //print_matrix(matrix, rows, cols);

    size_t res_rows, res_cols;
    double** result = decomposite_THEM(matrix, cols, rows, n, m, size, rank, &res_rows, &res_cols);
    //printf("\n");
    //print_matrix(result, res_rows, res_cols);

    free_matrix(&matrix, rows, cols);
    free_matrix(&result, res_rows, res_cols);
    return true;
}

bool test_transfer_data(int rows, int cols,int rank, int size)
{

    double** matrix = allocate_matrix(rows, cols);
    int i, j;
    for(i = 0; i < rows; ++i)
        for(j = 0; j < cols; ++j)
            matrix[i][j] = rank;

    double** result = transfer_data(matrix, rows, cols, rank, size);

    printf("%d\n", rank);
    print_matrix(result, rows, cols);

    return true;
}

bool is_transposed(double** m1, double** m2, size_t rows, size_t cols)
{
    int i, j;
    for(i = 0; i < rows; ++i)
        for(j = 0; j < cols; ++j)
            if(m1[i][j] != m2[j][i])
                return false;
    return true;
}

bool test_transpose(int rows, int cols, int rank, int size)
{
    int i, j;
    double** matrix = NULL;
    if(rank == 0)
    {
        //matrix = give_temp_matrix(rows, columns);
        matrix = give_temp_matrix(rows, cols);
      //  print_matrix(matrix, cols, rows);
    }
    else
    {
        //receive rows and columns
    }
    size_t new_rows;
    size_t new_cols = rows;
    double** result = DistributeColumns(matrix, cols, rows, size, rank, 0, &new_rows);

    //print_matrix(result, new_rows, rows);
      //  free_matrix(&matrix, columns, rows);
    //if(rank == 0)
    //  free_matrix(&matrix, cols, rows);//  print_matrix(matrix, columns, rows);
    //printf("-%d Entering transpose\n", rank);
    Transpose(&result, &new_cols, &new_rows, cols, rows, size, rank);
    //printf("--- %d , %d ===\n", rows, new_rows);
    //double** collected = allocate_matrix(new_rows, rows);
    //matrix = allocate_matrix(rows, cols);
    //printf("-%d left transpose\n", rank);
    //printf("%d \n", rank);
    //print_matrix(result, new_rows, rows);
//
    double** transposed = NULL;
    if(rank == 0)
        transposed = allocate_matrix(rows, cols);
    CollectColumns(result, new_rows, new_cols, size, rank, 0, transposed, rows, cols);

//    printf("\n");
//    if(rank == 0){ printf("%d %d\n", rows, cols);
//        print_matrix(matrix, rows, cols);}
    //if(rank == 0)
    //    print_matrix(transposed, rows, cols);

    if(rank == 0 && is_transposed(matrix, transposed, cols, rows))
        printf("%d, %d - ok\n", rows, cols);

    free_matrix(&result, new_rows, cols);

    if(rank == 0)
    {
        free_matrix(&matrix, cols ,rows);
        free_matrix(&transposed, rows, cols);
    }
}

int main (int argc, char* argv[])
{
    int rank, size, ierr, i, j, starter=0;
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    //test_distribute(3, 4, rank, size);
    //test_collect(101, 420, rank, size);
    //test_decomposite_me(34, 7, 100, 7, rank, size);
    //test_decomposite_them(3, 7, 6, 7, rank, size);
    //test_transfer_data(size, 1, rank, size);
    for(i = 0 ; i < 1000; ++i)
        for(j = 0 ; j < 1000; ++j)
            test_transpose(i, j, rank, size);

            //test_transpose(i, j, rank, size);
//    double** original_matrix = NULL;
//    size_t m = 4;
//    size_t n = 7;
//
//    if(rank == 0)
//    {
//        // the original matrix is transposed to pass easier
//        original_matrix = allocate_matrix(n, m);
//        int i, j;
//        for(i = 0; i < m; ++i)
//            for(j = 0; j < n; ++j)
//                original_matrix[j][i] = i*n + j;
//        // read matrix
//        // send matrix size
//        print_matrix(original_matrix, n, m);
//    }
//    else
//    {
//        // receive matrix size
//    }
//    //printf("%d %d " , rank , size);
//    //matrix = allocate_matrix(m, (int)(n/size - 0.5) + 1);
//    //print_matrix(matrix, m, (int)(n/size - 0.5) + 1);
//    size_t rows = -1, columns = m;
//    double** matrix = NULL;
//    matrix = DistributeColumns(original_matrix, m, n, size, rank, 0, &rows);
//
//    if(matrix == NULL)
//        exit(0);
//
//    size_t t_rows = rows, t_columns = m;
//    Transpose(&matrix, &t_columns, &t_rows, n, m, size, rank);
//
//    //print_matrix(matrix, t_rows, t_columns);
//    double** transpoed = allocate_matrix(m, n);
//    CollectColumns(transpoed, n, m, size, rank, 0, matrix);
//    if(rank == 0)
//    {
//        print_matrix(transpoed, m, n);
//        free_matrix(&original_matrix, n, m);
//        free_matrix(&transpoed, m, n);
//    }
//    free_matrix(&matrix, t_rows, t_columns);
//
//    matrix = allocate_matrix(4, 4);
//    i = 0;
//    for(; i < 4; ++i)
//    {
//        matrix[0][i] = rank*5 + i;
//        matrix[1][i] = rank*5 + i;
//        matrix[2][i] = rank*5 + i;
//        matrix[3][i] = rank*5 + i;
//    }
//    for(i =0 ; i < 4; ++i)
//    {
//        //printf("--%d, %d, %d, %d, %d \n", rank, (int)matrix[i][0], (int)matrix[i][1], (int)matrix[i][2], (int)matrix[i][3]);
//    }
//    double** res = transfer_data(matrix, 4, 4, rank, size);
//
//    for(i =0 ; i < 4; ++i)
//    {
//        printf("%d, %d, %d, %d, %d \n", rank, (int)res[i][0], (int)res[i][1], (int)res[i][2], (int)res[i][3]);
//    }
//

    MPI_Finalize();

    return 0;
}
