#ifndef SSOR
#define SSOR
#include "matrix_allocation.h"
#include "/usr/include/mpi/mpi.h"
//#include "mpi/mpi.h"

enum SendRecvPosition
{
    SEND_TOP = 0,
    RECV_TOP = 1,
    SEND_BOTTOM = 2,
    RECV_BOTTOM = 3,
    SEND_RIGHT = 4,
    RECV_RIGHT = 5,
    SEND_LEFT = 6,
    RECV_LEFT = 7
};

struct SendRecv
{
    int positions[8];
};

struct Info
{
    int N, P, Q, f; // f not really here /// PxQ(rows x cols)
    double w;
    double h;
    int l, m;
    int coords[2];
    int rank, size;
    MPI_Comm cartcomm;
    int first_cell; // 0 red, 1 black
    struct SendRecv data;
    double epsilon;
};

enum Color{ RED=0, BLACK=1};

struct SendRecvDatatypes
{
    MPI_Datatype stride_2;
    MPI_Datatype stride_2m;
};

struct Timing
{
    double time[5];
};

double g(double x, double y, double def);

void set_up_data_position(struct Matrix1D* u, struct Info* info);

void initialize_starting_values(struct Matrix1D* u, struct Info* info);

int is_cell_black(size_t i, size_t j);

struct SendRecvDatatypes createDataTypes(struct Matrix1D* u);
// returns the original most upper left coordinate of the info processor
// the returened coordinates are in x and y
void get_original_starting_coordinates(struct Info* info, int* x, int *y);

int absolute(int x);

//void create_datatype(int first_cell, int is_black, struct Datatype* data, int amount, MPI_Datatype* old);
//struct Transfer_data create_sending_data(struct Matrix1D* u, struct Info* info, enum Color color);

void transfer(struct Matrix1D* u, struct Info* info, enum Color color);

void set_up_grid(struct Info* info);

void solve(struct Matrix1D* u, struct Matrix1D* b, struct Info *info, struct Timing* timings, int max_steps);
void evaluate(struct Matrix1D* u, struct Info* info, struct Matrix1D* b, int is_black);

void calculate(struct Matrix1D* u, size_t i, size_t j, struct Info* info, struct Matrix1D* b); // correct b

double sum_up_all(struct Matrix1D* u);
double get_euclidian(int collector, struct Matrix1D* u);

#endif
