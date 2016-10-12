#ifndef CG
#define CG
#include "matrix_allocation.h"
#include "SSOR.h"
#include "/usr/include/mpi/mpi.h"

struct CG_Data
{
    struct Matrix1D x;
    struct Matrix1D s;
    struct Matrix1D r;
    struct Matrix1D As; // q
    double alpha;
    double beta;
    double ro;
};

/*
struct CG_Info
{
    int N, P, Q, f; // f not really here /// PxQ(rows x cols)
    double h;
    int l, m;
    int coords[2];
    int rank, size;
    MPI_Comm cartcomm;
    int first_cell; // 0 red, 1 black
    struct SendRecv data;
}; */

void AXPY(double alpha, struct Matrix1D* x, struct Matrix1D* y, struct Matrix1D* result, struct Info* info);

void XPb(struct Matrix1D* x, double b, struct Matrix1D* result, struct Info* info);

double DOT(struct Matrix1D* x, struct Matrix1D* y, struct Info* info);

void laplace_matrix_vector_product(struct Matrix1D* q, struct Matrix1D* s,struct Info *info);

void cg_evaluate(struct Matrix1D* u, struct Info* info, double b, int is_black);

void cg_calculate(struct Matrix1D* u, size_t i, size_t j, struct Info* info, double b);

void CG_solver(struct CG_Data* data, struct Info* info);

void allocate(struct CG_Data* data, struct Info* info);

void cg_settup();


#endif
