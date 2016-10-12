#ifndef MATRIX_ALLOCATION
#define MATRIX_ALLOCATION

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

struct Matrix1D
{
    double* data;
    size_t n; // rows
    size_t m; // columns
};

void copy_Matrix1D(struct Matrix1D* mat, struct Matrix1D* source);

bool allocate_Matrix1D(size_t rows, size_t columns, struct Matrix1D* matrix);

void free_Matrix1D(struct Matrix1D* matrix);

double* get_Matrix1D(struct Matrix1D* matrix, size_t i, size_t j);

void print_Matrix1D(struct Matrix1D* matrix);

double** allocate_matrix(size_t rows, size_t columns);

void free_matrix(double*** matrix, size_t rows, size_t columns);

void print_matrix(double** matrix, size_t rows, size_t columns);

void print_vector(double* vec, size_t elements);

double* matrix_to_vector(double** matrix, size_t rows, size_t columns);

double** vector_to_matrix(double* arr, size_t rows, size_t columns);

#endif
