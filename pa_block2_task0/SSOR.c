#include "SSOR.h"

//int g(size_t x, size_t y, int def)
//{
//    if(y = 0)
//        return 1 - x;
//    if(y = 1)
//        return x;
//    if(x = 0)
//        return 1 - y;
//    if(x = 1)
//        return y;
//    return def;
//}
//
//void initialize_starting_values(struct Matrix1D* u)
//{
//    memset(u->data, 0, sizeof(double));
//
//    // set the edges to g(x, y);
//
//}
//
//int is_cell_black(size_t i, size_t j)
//{
//    return ((i + j) % 2);
//}
//
//int absolute(int x)
//{
//    return (x > 0) ? x : -x;
//}
//
//void create_datatype(int first_cell, int is_black, struct Datatype* data, int amount, MPI_Datatype* old)
//{
//    if((first_cell == 0 && is_black == 1) || (first_cell == 1 && is_black == 0)) // abs(first_cell - is_black) == 1
//    {
//        data->start += amount;
//        data->len--;
//    }
//    MPI_Type_vector((int)(data->len / 2.0 + 0.5), 1, amount, *old, &data->type);
//    MPI_Type_commit(&data->type);
//}
//
//struct Transfer_data create_sending_data(struct Matrix1D* u, struct Info* info, enum Color color)
//{
//    struct Transfer_data td;
//
//    MPI_Datatype old;
//    MPI_Type_contiguous(1, MPI_DOUBLE, &old);
//    MPI_Type_commit(&old);
////    MPI_Type_vector(2, 1, 3, type2, &type);
////    MPI_Type_commit(&type);
//    int is_black = (color == BLACK) ? 1 : 0; // or just change everu is_black to color !!
//
//    int first_cell = info->first_cell;
//    td.top.len = u->m;
//    td.top.start = 0;
//    create_datatype(first_cell, is_black, &td.top, 1, &old);
//
//    first_cell = (u->n % 2 == 0) ? absolute(first_cell - 1) : first_cell;
//    td.bottom.len = u->m;
//    td.bottom.start = u->m * (u->n - 1); // u->n > 0 !!!
//    create_datatype(first_cell, is_black, &td.bottom, 1, &old);
//
//    first_cell = info->first_cell;
//    td.left.len = u->n;
//    td.left.start = 0;
//    create_datatype(first_cell, is_black, &td.left, u->m, &old);
//
//    first_cell = (u->n % 2 == 0) ? absolute(first_cell - 1) : first_cell;
//    td.right.len = u->n;
//    td.right.start = u->m - 1;
//    create_datatype(first_cell, is_black, &td.right, u->m, &old);
//
//    return td;
//}
//
//void transfer(struct Matrix1D* u, struct Info* info, enum Color color)
//{
//    int bottom, top, right, left;
//    struct  Transfer_data td = create_sending_data(u, info, color);
//    MPI_Status status;
//
//    MPI_Cart_shift(info->cartcomm, 0, 1, &top, &bottom);
//    MPI_Cart_shift(info->cartcomm, 0, 1, &left, &right);
//
//    // exchange rowwise
//    if(info->coords[0] % 2 == 0)
//    {
//        if(bottom > -1)
//        {
//            MPI_Send(u + td.bottom.start, 1, td.bottom.type, bottom, 0, MPI_COMM_WORLD); // maybe cartcomm
//            MPI_Recv(u + td.bottom.start, 1, td.bottom.type, bottom, 0, MPI_COMM_WORLD, &status);
//        }
//        if(top > -1)
//        {
//            MPI_Recv(u + td.top.start, 1, td.top.type, top, 0, MPI_COMM_WORLD, &status);
//            MPI_Send(u + td.top.start, 1, td.top.type, top, 0, MPI_COMM_WORLD); // maybe cartcomm
//        }
//    }
//    // missing else
//
//    // exchange columnwise
//    if(info->coords[1] % 2 == 0)
//    {
//        if(left > -1)
//        {
//            MPI_Send(u + td.left.start, 1, td.left.type, left, 0, MPI_COMM_WORLD); // maybe cartcomm
//            MPI_Recv(u + td.left.start, 1, td.left.type, left, 0, MPI_COMM_WORLD, &status);
//        }
//        if(right > -1)
//        {
//            MPI_Recv(u + td.right.start, 1, td.right.type, right, 0, MPI_COMM_WORLD, &status);
//            MPI_Send(u + td.right.start, 1, td.right.type, right, 0, MPI_COMM_WORLD); // maybe cartcomm
//        }
//    }
//    /// missing else
//}
//
//
//void set_up_grid(struct Info* info)
//{
//    int error;
//    enum N { ndims = 2 }; // zwi dimension
//    int dims[ndims] = {info->P , info->Q}; // {P, Q}
//    int periods[ndims] = {0, 0}; // this is fine no wrapping
//    int reorder = 0; // 0 for now
//    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &info->cartcomm);
//
//    error = MPI_Dims_create(info->size, ndims, dims);
//
//    error = MPI_Cart_coords(info->cartcomm, info->rank, ndims, info->coords);
//}
//
//void solve(struct Matrix1D* u,struct Info *info)
//{
//    int i, j=0;
//    bool critirion_off = true;
//    info->first_cell = is_cell_black(info->coords[0]  * info->P, info->coords[1] * info->Q);
//    while(critirion_off)
//    {
//        evaluate(u, info, 0, 0); // b is needed red
//
//        transfer(u, info, RED);
//
//        evaluate(u, info, 0, 1); // b is needed black
//        evaluate(u, info, 0, 1); // b is needed black
//
//        transfer(u, info, BLACK);
//
//        evaluate(u, info, 0, 0); // b is needed red
//
//        if(false) // stopping criterion
//            critirion_off = false;
//    }
//
//}
//
//void evaluate(struct Matrix1D* u, struct Info* info, double b, int is_black) // correct b
//{
//    int i, j, i_bound = u->n - 1, j_bound = u->m - 1, additional = is_black;
//
//    for(i = 1; i < i_bound; ++i)
//        for(j = 1 + additional - info->first_cell; j < j_bound; j += 2) // not going to work - firstcell
//            calculate(u, i, j, info, b); // pass b
//}
//
//void calculate(struct Matrix1D* u, size_t i, size_t j, struct Info* info, double b) // correct b
//{
//    double uij = *get_Matrix1D(u, i, j);
//    // the bounderies are set in the initializing step so don't worry about it.
//    uij = (info->w/4.0) * (info->h*info->h * b + *get_Matrix1D(u, i, j - 1) +
//                                                 *get_Matrix1D(u, i, j + 1) +
//                                                 *get_Matrix1D(u, i - 1, j) +
//                                                 *get_Matrix1D(u, i + 1, j)) + (1 - info->w) * uij;
//}
