#ifndef GAUSSIAN
#define GAUSSIAN

#include "matematica.h"
void interchange_rows(double** self, size_t rows, size_t columns, size_t first, size_t second)
{
    //printf("interchange row %d with row %d \n", first, second);
    if(first == second)
        return;
    size_t i;
    for(i = 0; i < rows; ++i)
        swap_double(&self[i][first], &self[i][second]);
}

size_t get_pivot_row(double* self, size_t length, size_t from)
{
    int current = from;
    size_t i;
    for(i = from + 1; i < length; ++i){
        if(absolute(self[i]) > absolute(self[current])){
            current = i;}}
    return current;
}

struct Package
{
    char* buff;
    int size;
    int position;
} ;

struct Matrix
{
    double** matrix;
    size_t rows;
    size_t cols;
};

void gaussian_elemination(double** self, size_t rows, size_t columns, size_t original_rows,
                             size_t original_columns, int rank, int size, int* P)
{
    size_t k = 0, bound, bound_k, i, j, e, data_size, buff_size;
    int position, pivot;
    double* data = NULL;
    char* buff = NULL;
    struct Package package;

    if(rank == 0)
        for(i = 0; i < original_columns; ++i)
            P[i] = i;

    bound_k = original_columns - 1;
    bound = bound_k + 1;
    for(k = 0; k < bound_k; ++k)
    {
        data_size = (columns - k - 1);
        data = malloc(data_size*sizeof(double) + sizeof(double));

        package.size = data_size * sizeof(double) + sizeof(int);
        package.buff = malloc(package.size);
        package.position = 0;

        if( k % size == rank)
        {
            e = (int)(k / size);
            // get pivot element
            pivot = get_pivot_row(self[e], columns ,k);
            // interchange rows
            interchange_rows(self, rows, columns, k, pivot);
            if(pivot == k && absolute(self[e][k]) <= 10e-7)
                memset(data, 0, sizeof(double) * (columns - k + 1));
            else
            {// calculate multipliers
                for(i = k + 1; i < bound; ++i)
                    self[e][i] = self[e][i] / self[e][k];

                for(j = 0, i = k + 1; i < bound; ++i, ++j) // memcpy
                    data[j] = self[e][i];
            }
            MPI_Pack(data, data_size, MPI_DOUBLE, package.buff, package.size, &package.position, MPI_COMM_WORLD);
            MPI_Pack(&pivot, 1, MPI_INT, package.buff, package.size, &package.position, MPI_COMM_WORLD);
        }
        MPI_Bcast(package.buff, package.size, MPI_PACKED, k % size, MPI_COMM_WORLD);

        if(k % size != rank)
        {
            MPI_Unpack(package.buff, package.size, &package.position, data, data_size, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Unpack(package.buff, package.size, &package.position, &pivot, 1, MPI_INT, MPI_COMM_WORLD);

            interchange_rows(self, rows, columns, k, pivot);
        }
        if(rank == 0 && pivot != k)
            swap_int(&P[k], &P[pivot]);

        e = my_rows_from(k + 1, rank, size);
        for(i = e; i < rows; ++i)
            for(j = k + 1; j < columns; ++j)
                self[i][j] = self[i][j] - data[j - k - 1] * self[i][k];

        free(data);
        free(package.buff);
    }
}

/// AX = B Anxn Bnxq solve for X vector with length n
double** triangular_solve(double** A, double** B, size_t n, size_t a_rows, size_t b_rows, size_t q,
                          int rank, int size)
{
    int to_rows = a_rows - 1;
    int current = (n-1)%size;
    int i = 0,j = 0, k = 0, e, bound, step = (q/size) + 1;
    double** result = allocate_matrix(a_rows, q);
    double* sums = malloc(sizeof(double) * (q + size));
    double* xes = malloc(sizeof(double) * q);
    double akk=1;

    for(i = n-1 ; i >= 0; --i)
    {
        e = i / size ;
        bound = e + ((rank <= current) ? 1 : 0); // my cols to
        for(k = 0; k < q; ++k)
            sums[k] = 0;

        for(j = a_rows - 1; j >= bound ; --j)
            for(k = 0; k < q; ++k)
                sums[k] += result[j][k] * A[j][i];

        //MPI_Reduce(sums, result[e], q, MPI_DOUBLE, MPI_SUM, current, MPI_COMM_WORLD);

        MPI_Allreduce(sums, xes, q, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /// fix this ///
        if(rank == current)
            akk = A[e][i];
        MPI_Bcast(&akk, 1, MPI_DOUBLE, current, MPI_COMM_WORLD);

        int at_x;
        for(at_x = rank, k = 0; k < b_rows; ++k, at_x += size)
            xes[k] = (B[k][i] - xes[at_x])/akk;

        MPI_Gather(xes, step, MPI_DOUBLE, sums, step, MPI_DOUBLE, current, MPI_COMM_WORLD);

        // this is kinda bad ..
        if(rank == current)
            for(k=0, at_x =0, j = 0; j < q; ++j)
            {
                result[e][j] = sums[at_x];
                at_x += step;
                if(at_x >= step*size)
                    at_x = ++k;
            }

        current = (current - 1 + size) % size;
    }
    free(xes);
    free(sums);
    return result;
}

/// AX = B Anxn Bnxq solve for X vector with length n
double** triangular_solve11(double** A, double** B, size_t n, size_t a_rows, size_t q,
                          int rank, int size)
{
    int to_rows = a_rows - 1;
    int current = (n-1)%size;
    int i = 0,j = 0, k = 0, e, bound;
    double** result = allocate_matrix(a_rows, q);
    double* sums = malloc(sizeof(double) * q);
    double* xes = malloc(sizeof(double) * (q/size));
    double akk;
    for(i = n-1 ; i >= 0; --i)
    {
        e = i / size ;
        bound = e + ((rank <= current) ? 1 : 0); // my cols to
        for(k = 0; k < q; ++k)
            sums[k] = 0;

        for(j = a_rows - 1; j >= bound ; --j)
            for(k = 0; k < q; ++k)
                sums[k] += result[j][k] * A[j][i];

        MPI_Reduce(sums, result[e], q, MPI_DOUBLE, MPI_SUM, current, MPI_COMM_WORLD);

        /// job done in one processor
        if(rank == current)
            for(k = 0; k < q; ++k)
                result[e][k] = (B[e][k] - result[e][k])/A[e][i];

        current = (current - 1 + size) % size;
    }
    free(xes);
    free(sums);
    return result;
}

double LUPTA_norm(double** LU, double** A, double** P, int n)
{
    double result = 0;
    double* AV = matrix_to_vector(A, n, n);
    double* L = matrix_to_vector(LU, n, n);
    double* PV = matrix_to_vector(P, n, n);

    int i,j;
    double* U = matrix_to_vector(LU, n, n);
    int bound = 0;
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < bound; ++j)
            U[i*n + j] = 0;
        ++bound;
    }

    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n , n, 1.0, L, n, U, n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1.0, PV, n, AV, n, 1.0, U, n);

    result = frobenius_norm(vector_to_matrix(U, n, n), n, n);
    free(L);
    free(AV);
    free(PV);
    free(U);
    return result;
}


double BUX_norm(double** B, double** U, double** X, int n, int q)
{
    double result = 0;

    double* BV = matrix_to_vector(B, n, q);
    double* UV = matrix_to_vector(U, n, n);
    double* XV = matrix_to_vector(X, n, q);

    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, q, 1.0, UV, n, XV, q);

    double* I = identity(n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, q, n, 1.0, I, n, BV, q, -1.0, XV, q);

    result = frobenius_norm(vector_to_matrix(XV, n, q), n, q);

    free(BV);
    free(UV);
    free(XV);
    free(I);
    return result;
}


#endif
