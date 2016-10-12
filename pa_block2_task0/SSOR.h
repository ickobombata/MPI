#ifndef SSOR
#define SSOR
#include <stdio.h>
//#include "matrix_allocation.h"
#include "/usr/include/mpi/mpi.h"



struct Info
{
    int N, P, Q, f; // f not really here
    double w;
    double h;
    int l, m;
    int coords[2];
    int rank, size;
    MPI_Comm cartcomm;
    int first_cell; // 0 red, 1 black
};

enum Color{ RED=0, BLACK=1};

struct Datatype
{
    MPI_Datatype type;
    size_t start;
    size_t len;
};

struct Transfer_data
{
    //double* left, *right, *top, *bottom;
    //size_t lsize, rsize, tsize, bsize;
    struct Datatype left;
    struct Datatype right;
    struct Datatype top;
    struct Datatype bottom;
};


int g(size_t x, size_t y, int def);

void initialize_starting_values(struct Matrix1D* u);

int is_cell_black(size_t i, size_t j);

//
//int absolute(int x);
//
//void create_datatype(int first_cell, int is_black, struct Datatype* data, int amount, MPI_Datatype* old);
//struct Transfer_data create_sending_data(struct Matrix1D* u, struct Info* info, enum Color color);
//
//void transfer(struct Matrix1D* u, struct Info* info, enum Color color);
//
//void set_up_grid(struct Info* info);
//
//void solve(struct Matrix1D* u,struct Info *info);
//void evaluate(struct Matrix1D* u, struct Info* info, double b, int is_black);
//
//void calculate(struct Matrix1D* u, size_t i, size_t j, struct Info* info, double b); // correct b


#endif
