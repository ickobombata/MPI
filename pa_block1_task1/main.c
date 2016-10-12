#include "/usr/include/mpi/mpi.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <cblas.h>
#include <unistd.h>

#include "matrix_operationMPI.h"
#include "gaussianMPI.h"
#define PI 3.14159265359

/*              m
matrix = [[1, 4, 7, 10],        1,   2,  3,
       n  [2, 5, 8, 11],        4,   5,  6,
          [3, 6, 9, 12]]   ==   7,   8,  9,
                                10, 11, 12.
*/

double** give_me_matrix_A(size_t n, size_t m, size_t* rows, size_t* columns, int rank, int size)
{
    *rows = my_rows_from(n, rank, size);
    *columns = n;
    double** result = allocate_matrix(*rows, *columns);
    size_t i,j;
    size_t new_j = rank;
    for(new_j = rank, j = 0; j < *rows; ++j, new_j += size)
    {
        for(i = 0; i < n ; ++i)
        {
            result[j][i] = ((i+new_j+2) / (n)) * sin(((i+1)*(new_j+1)*PI) / (n+1));
        }
     }
    return result;
}

double** give_me_matrix_B(size_t q, size_t n, size_t* rows, size_t* columns, int rank, int size)
{
    *rows = my_rows_from(q, rank, size);
    *columns = n;
    double** result = allocate_matrix(*rows, *columns);
    size_t i,j, mul = 0;
    size_t new_j;
    for(new_j = rank, j = 0; j < *rows; ++j, new_j += size)
    {
        for(i = 0; i < n ; ++i)
        {
            result[j][i] = i+new_j;
        }
     }
    return result;
}

void test(int n, int q, double* a, double* b, double* c)
{
    int rank, size, ierr, i, j, starter=0;
    MPI_Status status;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    double t1,t2 ,t3,t4, average1=0, average2=0;
    if(rank == 0)
        t1 = MPI_Wtime();

    double** X;
    double** matrix;

    size_t a_rows, b_rows, a_columns, b_columns;

    double** A = give_me_matrix_A(n, n, &a_rows, &a_columns, rank, size);
    double** B = give_me_matrix_B(q, n, &b_rows, &b_columns, rank, size);

    size_t x_rows = a_rows, x_cols = q;

    int* Per = NULL;

    if(rank == 0)
        Per = malloc(sizeof(int) * n);

    concatenate_matrix(&A, a_rows, n, B, b_rows);

    gaussian_elemination(A, a_rows + b_rows, n, n+q, n , rank, size, Per);

    if(rank == 0)
        t2 = MPI_Wtime();

    X = triangular_solve(A, B, n, a_rows, b_rows, q, rank ,size);

    if(rank == 0)
    {
        t3 = MPI_Wtime();
        *a = t2-t1;
        *b = t3-t1;
    }
    Transpose(&B, &b_columns, &b_rows, q, n, size, rank);

    X = triangular_solve11(A, B, n, a_rows,  q, rank ,size);

    if(rank == 0)
    {
        t4 = MPI_Wtime();
        *c = t4-t3 + t2-t1;
    }
    free_matrix(&A, a_rows, a_columns);
    free_matrix(&B, b_rows, b_columns);
}



int main (int argc, char* argv[])
{
    int rank, size, ierr, i, j, starter=0;
    MPI_Status status;

    int times = 5;
    int n = 500, q = 1000;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    double gauss, triangular_for_q, triangular_for_n;
    double sum[3] = {0.0, 0.0, 0.0};
    for(i=0; i < times; ++i)
    {
        test(n , q, &gauss, &triangular_for_q, &triangular_for_n);
        sum[0] += gauss; sum[1] += triangular_for_q; sum[2] += triangular_for_n;
    }
    if(rank == 0)
    {
    sum[0] /= times;
    sum[1] /= times;
    sum[2] /= times;
    printf("Time for gaussian :     %2.10lf\n\
            Time for q triangular : %2.10lf\n\
            Time for n triangular : %2.10lf\n", sum[0], sum[1], sum[2]);
    }
    MPI_Finalize();

    return 0;
}

//MPI_Comm_rank (MPI_COMM_WORLD, &rank);
//    MPI_Comm_size (MPI_COMM_WORLD, &size);
//    double** X;
//    double** matrix;
//
//    size_t a_rows, b_rows, a_columns, b_columns, q=1000, n=40;
//    size_t x_rows = a_rows, x_cols = q;
//    double** A = give_me_matrix_A(n, n, &a_rows, &a_columns, rank, size);
//    double** B = give_me_matrix_B(q, n, &b_rows, &b_columns, rank, size);
//    x_rows = a_rows;
//    size_t b0_rows = b_rows, b0_columns = b_columns;
//    size_t aO_rows = a_rows, aO_columns = a_columns;
//
//    double** A_original = give_me_matrix_A(n, n, &a_rows, &a_columns, rank, size);
//    double** B_original = give_me_matrix_B(q, n, &b_rows, &b_columns, rank, size);
//    int* Per = NULL;
//
//    if(rank == 0)
//        Per = malloc(sizeof(int) * n);
//
//    concatenate_matrix(&A, a_rows, n, B, b_rows);
//
//    gaussian_elemination(A, a_rows + b_rows, n, n+q, n , rank, size, Per);
//    double t1,t2 , average1=0, average2=0;
//
//    if(rank == 0)
//        t1 = MPI_Wtime();
//
//    X = triangular_solve(A, B, n, a_rows, b_rows, q, rank ,size);
//
//    if(rank == 0)
//        t2 = MPI_Wtime();
//
//    if(rank == 0)
//        printf("Time for first solver  %2.10lf\n", t2-t1);
//
//    if(rank == 0)
//        t1 = MPI_Wtime();
//
//    Transpose(&B, &b_columns, &b_rows, q, n, size, rank);
//
//    X = triangular_solve11(A, B, n, a_rows,  q, rank ,size);
//
//    if(rank == 0)
//        t2 = MPI_Wtime();
//
//    if(rank == 0)
//        printf("Time for second solver %2.10lf\n", t2-t1);
//
//    if(rank == 0)
//        matrix = allocate_matrix(n , q);
//
//    CollectColumns(X, a_rows, q, size, rank, 0, matrix, n, q);
//
//    Transpose(&B_original, &b0_columns, &b0_rows, q, n, size, rank);
//  //  Transpose(&B, &b_columns, &b_rows, q, n, size, rank);
//    Transpose(&A, &a_columns, &a_rows, n, n, size, rank);
//    Transpose(&A_original, &aO_columns, &aO_rows, n, n, size, rank);
//    //Transpose(&X, &x_cols, &x_rows, n, q, size, rank);
//
//    // collect LU
////    double** X0 = NULL;
////    double** B0 = NULL;
////    double** A0 = NULL;
////    double** Aoriginal0 = NULL;
////    double** Bvylna0 = NULL;
////
////    if(rank == 0)
////    {
////        X0 = allocate_matrix(n, q);
////        B0 = allocate_matrix(n, q);
////        A0 = allocate_matrix(n, n);
////        Aoriginal0 = allocate_matrix(n, n);
////        Bvylna0 = allocate_matrix(n, q);
////
////    }
////    CollectColumns(A_original, a_rows, a_columns, size, rank, 0, Aoriginal0, n, n);
////    CollectColumns(A, a_rows, a_columns, size, rank, 0, A0, n, n);
////    CollectColumns(B, b_rows, b_columns, size, rank, 0, Bvylna0, n, q);
////    CollectColumns(B_original, b0_rows, b0_columns, size, rank, 0, B0, n, q);
////    CollectColumns(X, x_rows, x_cols, size, rank, 0, X0, n, q);
////
////    if(rank == 0)
////    {
////        double** P = construct_permutation_matrix(Per, n);
////
////        /// HERE CODE FOR THE NORMS
////
////        //print_matrix(Aoriginal0, n, n);
////        //printf("\n");
////        //print_matrix(Bvylna0, n, q);
////
////        printf("LUPTA=%2.10lf, ", LUPTA_norm(A0, Aoriginal0, P, n));
////
////
////        printf("BUX=%2.10lf, \n", BUX_norm(Bvylna0, A0, X0, n, q));
////
////        free_matrix(&A0, n, n);
////        free_matrix(&Aoriginal0, n, n);
////        free_matrix(&Bvylna0, n, q);
////        free_matrix(&B0, n, q);
////        free_matrix(&X0, n, q);
////    }
//
//    if(rank == 0)
//    {
//        //print_matrix(matrix, n, q);
//        free_matrix(&matrix, n, q);
//    }
//
//    //Transpose(&A, &a_columns, &a_rows, n, n, size, rank);
//    free_matrix(&A, a_rows, a_columns);
//    free_matrix(&X, x_rows, q);
//    free_matrix(&B, b_rows, b_columns);
//



//    double** ca,** cb;
//    if(rank == 0)
//    {
//        ca = allocate_matrix(n, n);
//        cb = allocate_matrix(q, n);
//    }
//    //printf("%d\n", rank);
//    //print_matrix(A, a_rows, a_columns);
//    CollectColumns(A, a_rows, a_columns, size, rank, 0, ca, n, n);
//    CollectColumns(B, b_rows, b_columns, size, rank, 0, cb, q, n);
//
//    if(rank == 0)
//    {
//        print_matrix(ca, n, n);
//        free_matrix(&ca, n, n);
//        printf("\n");
//        print_matrix(cb, q, n);
//        free_matrix(&cb, q, n);
//        printf("\n");
//    }

