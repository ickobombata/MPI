#include <math.h>
#include <stdio.h>
#include "SSOR.h"

double g(double x, double y, double def)
{
    if(y == 0)
        return 1. - x;
    if(y == 1) //hmmm
        return x;
    if(x == 0)
        return 1. - y;
    if(x == 1)//hmmm
        return y;
    return def;
}

void initialize_starting_values(struct Matrix1D* u, struct Info* info)
{
    int i, bound_i = u->m * u->n;
    for(i = 0; i < bound_i ; ++i)
        u->data[i] = 0.;

    int x, y;
    int n_times_m = u->m * u->n;
    get_original_starting_coordinates(info, &x, &y);
    double bound;
    if(x == 0)
    {
        for(i = 0, bound = y*info->h ; i < n_times_m; i+= u->m, bound += info->h)
        {
            u->data[i] = g(0, bound, 0);
                //printf("calc %2.2lf %d %d\n", u->data[i], bound, i);
        }
    }
    if(x == info->N - info->m)
    {
        double new_x = 1;
        for(i = u->m - 1 ,bound = y*info->h; i < n_times_m; i += u->m, bound += info->h)
            u->data[i] = g(new_x, bound, 0);
    }

    if(y == 0)
    {
        for(i = 0 , bound = x*info->h; i < u->m; ++i, bound += info->h)
            u->data[i] = g(bound, 0, 0);
    }
    if(y == info->N - info->l)
    {
        double new_y = 1;
        for(i = u->m * u->n - u->m, bound = x*info->h ; i < n_times_m; ++i, bound += info->h)
        {
            u->data[i] = g(bound, new_y, 0);
            //printf("calc %2.2lf %d %d\n", u->data[i], bound, i);
        }
    }
}


void set_up_data_position(struct Matrix1D* u, struct Info* info)
{
    info->data.positions[SEND_TOP] = 1 + u->m;
    info->data.positions[RECV_TOP] = 1;

    info->data.positions[SEND_BOTTOM] = u->m * (u->n - 2) + 1;
    info->data.positions[RECV_BOTTOM] = u->m * (u->n - 1) + 1;

    info->data.positions[SEND_RIGHT] = u->m + u->m - 2;
    info->data.positions[RECV_RIGHT] = u->m + u->m - 1;

    info->data.positions[SEND_LEFT] = u->m + 1;
    info->data.positions[RECV_LEFT] = u->m;
}

void get_original_starting_coordinates(struct Info* info, int* x, int *y)
{
    *y = info->coords[0] * info->l;
    *x = info->coords[1] * info->m;
}

void set_up_grid(struct Info* info)
{
    int error;
    enum N { ndims = 2 }; // zwеi dimension
    int dims[ndims] = {info->P , info->Q}; // {P, Q}
    int periods[ndims] = {0, 0}; // this is fine no wrapping
    int reorder = 0; // 0 for now

    // creates a PxQ grid communication

    error = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &info->cartcomm);

    // distributes the processors on the grid
    error = MPI_Dims_create(info->size, ndims, dims);

    // creates cartesian coordinates for the processor on the grid
    error = MPI_Cart_coords(info->cartcomm, info->rank, ndims, info->coords);


    int x, y;
    get_original_starting_coordinates(info, &x, &y);
    info->first_cell = is_cell_black(x, y);
}

int is_cell_black(size_t i, size_t j)
{
    return ((i + j) % 2);
}

int absolute(int x)
{
    return (x > 0) ? x : -x;
}

struct SendRecvDatatypes createDataTypes(struct Matrix1D* u)
{
    struct SendRecvDatatypes srdt;

    MPI_Datatype old;
    MPI_Type_contiguous(1, MPI_DOUBLE, &old);
    MPI_Type_commit(&old);

    MPI_Type_vector((int)((u->m - 2) / 2.0 + 0.5), 1, 2, old, &srdt.stride_2);
    MPI_Type_commit(&srdt.stride_2);

    MPI_Type_vector((int)((u->n - 2) / 2.0 + 0.5), 1, 2*u->m, old, &srdt.stride_2m);
    MPI_Type_commit(&srdt.stride_2m);

    return srdt;
}

//////////////////////////////////
//////////////////////////////////
//
//int start_top_send(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    return 1 + u->m + ((absolute(first_cell - color)) ? 1 : 0);
//}
//int start_top_recv(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    return 1 + ((absolute(first_cell - color)) ? 0 : 1);
//}
//int start_bottom_send(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    int bottom_cell_color = ((u->n % 2) ? (first_cell) : (absolute(first_cell - 1)));
//    return u->m * (u->n - 2) + 1 + ((absolute(bottom_cell_color - color)) ? 1 : 0);
//}
//int start_bottom_recv(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    int bottom_cell_color = ((u->n % 2 == 0) ? (first_cell) : absolute(first_cell - 1));
//    return u->m * (u->n - 1) + 1 + ((absolute(bottom_cell_color - color)) ? 1 : 0);
//}
//
//int start_right_send(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    int right_cell_color = ((u->m % 2 == 0) ? absolute(first_cell - 1) : first_cell);
//    return u->m + u->m - 2 + ((absolute(right_cell_color - color)) ? u->m : 0);
//}
//int start_right_recv(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    int right_cell_color = ((u->m % 2 == 0) ? first_cell : absolute(first_cell - 1));
//    return u->m + u->m - 1 + ((absolute(right_cell_color - color)) ? u->m : 0);
//}
//
//int start_left_send(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    return u->m + 1 + ((absolute(first_cell - color) ? u->m : 0));
//}
//
//int start_left_recv(struct Matrix1D* u, enum Color color, int first_cell)
//{
//    return u->m + ((absolute(first_cell - color) ? 0 : u->m));
//}


double* top_send(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[SEND_TOP] + (color ^ info->first_cell);
}

double* top_recv(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[RECV_TOP] + !(color ^ info->first_cell);
}

double* bottom_send(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[SEND_BOTTOM] + (color ^ ((u->n & 1) ? info->first_cell : info->first_cell ^ 1));
}

double* bottom_recv(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[RECV_BOTTOM] + (color ^ ((u->n & 1) ? info->first_cell ^ 1 : info->first_cell));
}

double* right_send(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[SEND_RIGHT] + ((color ^ ((u->m & 1) ? info->first_cell : info->first_cell ^ 1)) ? u->m : 0);
}

double* right_recv(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[RECV_RIGHT] + ((color ^ ((u->m & 1) ? info->first_cell ^ 1 : info->first_cell)) ? u->m : 0);
}

double* left_send(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[SEND_LEFT] + ((color ^ info->first_cell) ? u->m : 0);
}

double* left_recv(struct Matrix1D* u, struct Info* info, enum Color color)
{
    return u->data + info->data.positions[RECV_LEFT] + ((color ^ info->first_cell) ? 0 : u->m);
}

///////////////////////////////
/////////////////////////
void transfer(struct Matrix1D* u, struct Info* info, enum Color color)
{
    int bottom, top, right, left;
    //struct  Transfer_data td = create_sending_data(u, info, color);
    MPI_Status status;

    MPI_Cart_shift(info->cartcomm, 0, 1, &top, &bottom);
    MPI_Cart_shift(info->cartcomm, 1, 1, &left, &right);

    struct SendRecvDatatypes srdt = createDataTypes(u);

    if(info->coords[0] % 2 == 0)
    {
        if(bottom > -1)
        {
            MPI_Send(bottom_send(u, info, color), 1, srdt.stride_2, bottom, 0, MPI_COMM_WORLD); // maybe cartcomm
            MPI_Recv(bottom_recv(u, info, color), 1, srdt.stride_2, bottom, 0, MPI_COMM_WORLD, &status);
        }
        if(top > -1)
        {
            MPI_Send(top_send(u, info, color), 1, srdt.stride_2, top, 0, MPI_COMM_WORLD); // maybe cartcomm
            MPI_Recv(top_recv(u, info, color), 1, srdt.stride_2, top, 0, MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        if(top > -1)
        {
            MPI_Recv(top_recv(u, info, color), 1, srdt.stride_2, top, 0, MPI_COMM_WORLD, &status);
            MPI_Send(top_send(u, info, color), 1, srdt.stride_2, top, 0, MPI_COMM_WORLD); // maybe cartcomm
        }
        if(bottom > -1)
        {
            MPI_Recv(bottom_recv(u, info, color), 1, srdt.stride_2, bottom, 0, MPI_COMM_WORLD, &status);
            MPI_Send(bottom_send(u, info, color), 1, srdt.stride_2, bottom, 0, MPI_COMM_WORLD); // maybe cartcomm

        }
    }
     //exchange columnwise
    if(info->coords[1] % 2 == 0)
    {
        if(left > -1)
        {
            MPI_Send(left_send(u, info, color), 1, srdt.stride_2m, left, 0, MPI_COMM_WORLD); // maybe cartcomm
            MPI_Recv(left_recv(u, info, color), 1, srdt.stride_2m, left, 0, MPI_COMM_WORLD, &status);
        }
        if(right > -1)
        {
            MPI_Send(right_send(u, info, color), 1, srdt.stride_2m, right, 0, MPI_COMM_WORLD); // maybe cartcomm
            MPI_Recv(right_recv(u, info, color), 1, srdt.stride_2m, right, 0, MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        if(right > -1)
        {
            MPI_Recv(right_recv(u, info, color), 1, srdt.stride_2m, right, 0, MPI_COMM_WORLD, &status);
            MPI_Send(right_send(u, info, color), 1, srdt.stride_2m, right, 0, MPI_COMM_WORLD); // maybe cartcomm
         }
        if(left > -1)
        {
            MPI_Recv(left_recv(u, info, color), 1, srdt.stride_2m, left, 0, MPI_COMM_WORLD, &status);
            MPI_Send(left_send(u, info, color), 1, srdt.stride_2m, left, 0, MPI_COMM_WORLD); // maybe cartcomm
        }
    }
}

void solve(struct Matrix1D* u, struct Matrix1D* b, struct Info *info, struct Timing* timings, int max_steps)
{
    int i, j=0;
    int critirion_off = 1, counter = 0;
    double norm[2] = {0.0, 0.0 };
    double epsilon = 0.000000000001;
    double sum_time[2] = {0.0, 0.0};
    // Default values.

    int normIdx = 0;
    int x, y;
    get_original_starting_coordinates(info, &x, &y);
    info->first_cell = is_cell_black(x, y);

   // for(i = 0 ; i < 1000; ++i)
    while(critirion_off && counter < max_steps)
    {
        evaluate(u, info, b, RED); // b is needed red

        transfer(u, info, RED);
        evaluate(u, info, b, BLACK); // b is needed blackprint_Matrix1D(u);
   // printf("\n");
        evaluate(u, info, b, BLACK); // b is needed black
//print_Matrix1D(u);
//printf("\n");
        transfer(u, info, BLACK);

        evaluate(u, info, b, RED); // b is needed red
//print_Matrix1D(u);
  //  printf("\n");


        MPI_Barrier(MPI_COMM_WORLD);
       // norm[normIdx] = sqrt(get_euclidian(0, u));

        MPI_Barrier(MPI_COMM_WORLD);
        sum_time[0] += MPI_Wtime();

        norm[normIdx] = sqrt(get_euclidian(0, u));

        if (info->rank == 0)
        {
            //printf("%10.10lf, %10.10lf\n", norm[0], norm[1]);
             if (fabs(norm[normIdx] - norm[normIdx ^ 1]) < epsilon)
             {
                critirion_off = 0;
             }
        }

        MPI_Bcast(&critirion_off, 1, MPI_INT, 0, MPI_COMM_WORLD);

        normIdx ^= 1;

        MPI_Barrier(MPI_COMM_WORLD);
        sum_time[1] += MPI_Wtime();
        ++counter;
    }
    timings->time[2] = sum_time[0]/counter;
    timings->time[3] = sum_time[1]/counter;
}

void evaluate(struct Matrix1D* u, struct Info* info, struct Matrix1D* b, int is_black) // correct b
{
    int i, j, i_bound = u->n - 1, j_bound = u->m - 1, additional = info->first_cell ^ is_black;

    for(i = 1; i < i_bound; ++i)
    {
        for(j = 1 + additional; j < j_bound; j += 2)
            calculate(u, i, j, info, b); // pass b
        additional = (additional ^ 1);
    }
}

void calculate(struct Matrix1D* u, size_t i, size_t j, struct Info* info, struct Matrix1D* b) // correct b
{
    double* uij = get_Matrix1D(u, i, j);

    // the bounderies are set in the initializing step so don't worry about it.
    *uij = (info->w/4.0) *
                ((double)info->h*info->h * *get_Matrix1D(b, i, j) +
                                               *get_Matrix1D(u, i, j - 1) +
                                               *get_Matrix1D(u, i, j + 1) +
                                               *get_Matrix1D(u, i - 1, j) +
                                               *get_Matrix1D(u, i + 1, j)) + (1.0 - info->w) * (*uij);
    //if(i == 1 && j == 1)
    //    printf("%2.2lf = %2.2lf, %2.2lf, %2.2lf, %2.2lf\n", *uij, *get_Matrix1D(u, i, j - 1), *get_Matrix1D(u, i, j + 1), *get_Matrix1D(u, i - 1, j) ,*get_Matrix1D(u, i + 1, j));
}

double get_euclidian(int collector, struct Matrix1D* u)
{
    double mine_sum = sum_up_all(u);
    double sum = 0;

    MPI_Reduce(&mine_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, collector, MPI_COMM_WORLD);

    return sum;
}

double sum_up_all(struct Matrix1D* u)
{
    int i,j;
    double sum = 0;
    for(i = 1; i < u->n - 1; ++i)
        for(j = 1; j < u->m - 1; ++j)
            sum += u->data[i * u->m + j] * u->data[i * u->m + j];
    return sum;
}
