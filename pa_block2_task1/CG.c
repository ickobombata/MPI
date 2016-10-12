#include "CG.h"

void AXPY(double alpha, struct Matrix1D* x, struct Matrix1D* y, struct Matrix1D* result, struct Info* info)
{
    int i, j, i_bound = x->n - 1, j_bound = x->m - 1;

    for(i = 1; i < i_bound; ++i)
    {
        for(j = 1; j < j_bound; ++j)
            *get_Matrix1D(result, i, j) = *get_Matrix1D(y, i, j) + alpha*(*get_Matrix1D(x, i, j));
    }
}

double DOT(struct Matrix1D* x, struct Matrix1D* y, struct Info* info)
{
    double result = 0.0, sum_all = 0.0;

    int i, j, i_bound = x->n - 1, j_bound = x->m - 1;

    for(i = 1; i < i_bound; ++i)
    {
        for(j = 1; j < j_bound; ++j)
            result += (*get_Matrix1D(x, i, j)) * (*get_Matrix1D(y, i, j));
    }

    MPI_Allreduce(&result, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sum_all;
}

void update(struct Matrix1D* q, struct Matrix1D* s,struct Info *info)
{
    int i, j;

    for(i = 1 ; i < q->n - 1 ; ++i)
        for(j = 1; j < q->m - 1; ++j)
            {*get_Matrix1D(q, i, j) = 4 * (*get_Matrix1D(s, i, j)) - *get_Matrix1D(s, i, j - 1)
                                                                  - *get_Matrix1D(s, i, j + 1)
                                                                  - *get_Matrix1D(s, i - 1, j)
                                                                  - *get_Matrix1D(s, i + 1, j);
            //if(info->rank == 1)
            //    printf("%d, %d, %2.2lf\n", i, j, *get_Matrix1D(q, i, j));
            }
}

void laplace_matrix_vector_product(struct Matrix1D* u, struct Matrix1D* s, struct Info *info)
{
    transfer(s, info, RED);
    transfer(s, info, BLACK);

//    if(info->rank == 1)
//    {
//        print_Matrix1D(s);
//        printf("\n");
//        print_Matrix1D(u);
//        printf("\n");
//    }
    update(u, s, info);
//    if(info->rank == 1)
//    {
//        print_Matrix1D(u);
//        printf("\n");
//    }
}
//
void zero_up(struct Matrix1D* u)
{
    int um = u->m * u->n, i = 0;
    for(; i < um; ++i)
        u->data[i] = 0;
}

void initialize_borders(struct Matrix1D* u, struct Info* info)
{
    int i, x, y;
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

void minus_matrix(struct Matrix1D* u)
{
    int i, j, i_bound = u->n - 1, j_bound = u->m - 1;

    for(i = 1; i < i_bound; ++i)
        for( j = 1; j < j_bound; ++j)
            *get_Matrix1D(u, i, j) = -*get_Matrix1D(u, i, j);
}
void set_borders(struct Matrix1D* u, struct Info* info, double val)
{
    int i, x, y;
    int n_times_m = u->m * u->n;
    get_original_starting_coordinates(info, &x, &y);
    double bound;
    if(x == 0)
    {
        for(i = 0, bound = y*info->h ; i < n_times_m; i+= u->m, bound += info->h)
        {
            u->data[i] = val;
                //printf("calc %2.2lf %d %d\n", u->data[i], bound, i);
        }
    }
    if(x == info->N - info->m)
    {
        double new_x = 1;
        for(i = u->m - 1 ,bound = y*info->h; i < n_times_m; i += u->m, bound += info->h)
            u->data[i] = val;
    }

    if(y == 0)
    {
        for(i = 0 , bound = x*info->h; i < u->m; ++i, bound += info->h)
            u->data[i] = val;
    }
    if(y == info->N - info->l)
    {
        double new_y = 1;
        for(i = u->m * u->n - u->m, bound = x*info->h ; i < n_times_m; ++i, bound += info->h)
        {
            u->data[i] = val;
            //printf("calc %2.2lf %d %d\n", u->data[i], bound, i);
        }
    }
}

void precondisioned_CG_solver(struct CG_Data* data, struct Info* info, int steps)
{
    int i;
    double denominator = 1.0;
    double ro_prev;
    struct Matrix1D z;
    struct Timing t;

    allocate_Matrix1D(data->x.n, data->x.m, &z);
    //ran(&data->x);
//    printf("\n\n");
//    print_Matrix1D(&data->x);
//    printf("\n");

    //zero_up(&data->As);
    zero_up(&data->r);
    zero_up(&data->x);
    zero_up(&data->s);
    zero_up(&data->As);
    zero_up(&z);

    //print_Matrix1D(&data->s);


    initialize_borders(&data->x, info);

    /// r0 = b - Ax0
    copy_Matrix1D(&data->s, &data->x);
    laplace_matrix_vector_product(&data->r, &data->x, info);
    minus_matrix(&data->r);

    /// solve the system Tz0 = r0
    //initialize_borders(&z, info);
    solve(&z, &data->r, info, &t, steps);

    /// s0 := z0
    copy_Matrix1D(&data->s, &z);

    //set_borders(&data->s, info, 0.0);
    //initialize_borders(&data->s, info);
    //initialize_borders(&data->As, info);
    //initialize_borders(&data->r, info);

//    if(info->rank == 0)
//    {
        //printf("\nIteration -1\n0\n");
//        print_Matrix1D(&data->r);
//        printf("\n");
        //print_Matrix1D(&z);
        //printf("\n");
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 1)
//    {
        //print_Matrix1D(&data->s);
        //printf("\n");
//    }
    /// dot product : r0 * z0
    ro_prev = (DOT(&data->r, &z, info));

//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 0)
//    {
//        printf("ro_prev = %2.2lf\n", ro_prev);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 1)
//    {
//        printf("ro_prev = %2.2lf\n", ro_prev);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
    for(i = 0 ; i < 1000; ++i)
    {
        //printf("i = %d\n", i);
//        copy_Matrix1D(&data->As, &data->s);
        /// qk = A* sk .aka. As = A*s
        laplace_matrix_vector_product(&data->As, &data->s, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nIteration %d\n", i);
//
            //printf("\nAs %d\n", i);
            //print_Matrix1D(&data->As);
            //printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
            //print_Matrix1D(&data->s);
            //printf("\n");
//        }
        /// denom = sk * A * sk
        denominator = DOT(&data->s, &data->As, info);
        //printf("denominator = %2.2lf\n", denominator);
        //printf("denom = %2.5lf\n", denominator);
        if(denominator == 0)
        {
            printf("!!!!! %d \n", i );
            return; // do something with thits
        }
        /// alpha = rok / denom
        data->alpha = ro_prev / denominator;
        //printf("alpha = %2.2lf\n", data->alpha);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//            printf("%2.4lf / %2.2lf = %2.2lf\n", ro_prev, denominator, data->alpha);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//            printf("%2.4lf / %2.2lf = %2.2lf\n", ro_prev, denominator, data->alpha);
        /// xk+1 = xk + alpha*sk
        AXPY(data->alpha, &data->s, &data->x, &data->x, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nX %d\n", i);
//            print_Matrix1D(&data->x);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->x);
//            printf("\n");
//        }
        /// rk+1 = rk - alpha * As
        AXPY(-data->alpha, &data->As, &data->r, &data->r, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nR %d\n", i);
//            print_Matrix1D(&data->r);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->r);
//            printf("\n");
//        }
        /// solve the system Tzk+1 = rk+1
        //printf("Entring solve \n");
        zero_up(&z);
        solve(&z, &data->r, info, &t, steps);
        //copy_Matrix1D(&z, &data->r);
        //printf("EXIT SOlve\n");
        /// rok+1 = rok*zk+1
        data->ro = (DOT(&data->r, &z, info));

//        if(sqrt(data->ro) < info->epsilon)
//        {
//            printf("\n---%d steps --- %2.10lf, %2.10lf \n", i, sqrt(data->ro), data->ro);
//            return;
//        }

        /// beta = rok+1 / rok
        data->beta = - (data->ro / ro_prev);
//        printf("beta = %2.2lf\n", data->beta);
//
//        printf("ro_prev = %2.2lf\n", ro_prev);
        ro_prev = data->ro;
//        printf("ro = %2.2lf\n", ro_prev);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//            printf("beta = %2.4lf,  ro_prev=%2.2lf ro=%2.2lf\n", data->beta, ro_prev, data->ro);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//            printf("beta = %2.4lf,  ro_prev=%2.2lf ro=%2.2lf\n", data->beta, ro_prev, data->ro);
        /// sk+1 = zk+1 - beta*sk
        AXPY(-data->beta, &data->s, &z, &data->s, info);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nS %d\n", i);
//            print_Matrix1D(&data->s);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->s);
//            printf("\n");
//        }
    }
}

void CG_solver(struct CG_Data* data, struct Info* info)
{
    int i;
    double denominator = 1.0;
    double ro_prev;
    //ran(&data->x);
//    printf("\n\n");
//    print_Matrix1D(&data->x);
//    printf("\n");

    //zero_up(&data->As);
    zero_up(&data->r);
    zero_up(&data->x);
    zero_up(&data->s);
    zero_up(&data->As);
    initialize_borders(&data->x, info);

    copy_Matrix1D(&data->s, &data->x);
    laplace_matrix_vector_product(&data->r, &data->x, info);
    //XPb(&data->r, info->f, &data->r, info);
    minus_matrix(&data->r);

    copy_Matrix1D(&data->s, &data->r);

    //initialize_borders(&data->s, info);
    //initialize_borders(&data->As, info);
    //initialize_borders(&data->r, info);

//    if(info->rank == 0)
//    {
//        printf("\nIteration -1\n0\n");
//        print_Matrix1D(&data->s);
//        printf("\n");
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 1)
//    {
//        print_Matrix1D(&data->s);
//        printf("\n");
//    }

    ro_prev = (DOT(&data->r, &data->r, info));

//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 0)
//    {
//        printf("ro_prev = %2.2lf\n", ro_prev);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info->rank == 1)
//    {
//        printf("ro_prev = %2.2lf\n", ro_prev);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
    for(i = 0 ; true ; ++i)
    {
        //printf("i = %d\n", i);
//        copy_Matrix1D(&data->As, &data->s);
        laplace_matrix_vector_product(&data->As, &data->s, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nIteration %d\n", i);
//
//            printf("\nAs %d\n", i);
//            print_Matrix1D(&data->As);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->As);
//        }

        denominator = DOT(&data->s, &data->As, info);
        //printf("denominator = %2.2lf\n", denominator);
        //printf("denom = %2.5lf\n", denominator);
        if(denominator == 0)
        {
            printf("!!!!! %d \n", i );
            return; // do something with thits
        }
        data->alpha = ro_prev / denominator;
        //printf("alpha = %2.2lf\n", data->alpha);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//            printf("%2.4lf / %2.2lf = %2.2lf\n", ro_prev, denominator, data->alpha);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//            printf("%2.4lf / %2.2lf = %2.2lf\n", ro_prev, denominator, data->alpha);

        AXPY(data->alpha, &data->s, &data->x, &data->x, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nX %d\n", i);
//            print_Matrix1D(&data->x);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->x);
//            printf("\n");
//        }

        AXPY(-data->alpha, &data->As, &data->r, &data->r, info);

//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nR %d\n", i);
//            print_Matrix1D(&data->r);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->r);
//            printf("\n");
//        }

        data->ro = (DOT(&data->r, &data->r, info));

        if(sqrt(data->ro) < info->epsilon)
        {
            printf("steps !! %d\n\n", i);
           return;
        }

        data->beta = - (data->ro / ro_prev);
//        printf("beta = %2.2lf\n", data->beta);
//
//        printf("ro_prev = %2.2lf\n", ro_prev);
        ro_prev = data->ro;
//        printf("ro = %2.2lf\n", ro_prev);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//            printf("beta = %2.4lf,  ro_prev=%2.2lf ro=%2.2lf\n", data->beta, ro_prev, data->ro);
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//            printf("beta = %2.4lf,  ro_prev=%2.2lf ro=%2.2lf\n", data->beta, ro_prev, data->ro);

        AXPY(-data->beta, &data->s, &data->r, &data->s, info);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 0)
//        {
//            printf("\nS %d\n", i);
//            print_Matrix1D(&data->s);
//            printf("\n");
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(info->rank == 1)
//        {
//            print_Matrix1D(&data->s);
//            printf("\n");
//        }
    }
}

void allocate(struct CG_Data* data, struct Info* info)
{
    allocate_Matrix1D(info->l+2, info->m+2, &data->x);
    allocate_Matrix1D(info->l+2, info->m+2, &data->r);
    allocate_Matrix1D(info->l+2, info->m+2, &data->As);
    allocate_Matrix1D(info->l+2, info->m+2, &data->s);
}

int testCompare(struct Matrix1D * matrix, struct Info* info, double h)
{
    int i,j;
    double hSqr = h*h;
    double epsilon = 1.0e-7;
    int xx=0, yy=0;
    get_original_starting_coordinates(info, &xx, &yy);
    double x, y;

    for(i = 1, y =(double)(1+ yy) * h; i < matrix->n - 1 ; ++i, y += h)
    {
        for(j = 1, x=(double)(1+xx) * h; j < matrix->m - 1; ++j, x += h)
        {
            if (!(fabs(*get_Matrix1D(matrix, i, j) - (2. * x * y - x - y + 1.)) < epsilon))
            {
                printf("%d , (%d, %d) %10.9lf == %10.9lf\n", info->rank, i, j, *get_Matrix1D(matrix, i, j),
                                                                             (2. * x * y - x - y + 1.));
                return 0;
            }
        }
    }

    return 1;
}

void qq(struct Matrix1D* matrix, double h)
{
    int i,j;
    for(i = 0; i < matrix->n; ++i)
    {
        for(j = 0; j < matrix->m; ++j)
        {
            printf("%6.2lf, ", (2. * i * h * h * j - i*h - j*h + 1.));
        }
        printf("\n");
    }
}


void cg_settup()
{
    int error;
    struct Info info = {480, 1, 1, 0, 1.0};       //N=24, 48, 120, P = 2, Q = 3, f ≡ 0
    info.h = 1.0 / (double)(info.N + 1);
    info.l = info.N / info.P;
    info.m = info.N / info.Q;
    info.epsilon = 1.0e-9;

    MPI_Comm_rank (MPI_COMM_WORLD, &info.rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &info.size);        /* get number of processes */

    //printf("%d, %d, %d, %d, %d, %2.2lf, %2.2lf, %d, %d \n", info.rank, info.N, info.P, info.Q, info.f, info.w, info.h, info.l, info.m);
    struct Timing timings;
    MPI_Barrier(MPI_COMM_WORLD);
    timings.time[0] = MPI_Wtime();

    struct CG_Data data;
    allocate(&data, &info); /// TODO check memory

    set_up_grid(&info);

    //initialize_starting_values(&data.x, &info);

    //initialize_starting_values(&data.s, &info);

    set_up_data_position(&data.x, &info);

//    MPI_Barrier(MPI_COMM_WORLD);
//    timings.time[1] = MPI_Wtime();
    int sor_steps = 1;
    int steps = 0, max_steps = 1;

    timings.time[0] = 0;
    timings.time[4] = 0;
    for(steps = 0; steps < max_steps; ++steps)
    {
        timings.time[0] += MPI_Wtime();

        //precondisioned_CG_solver(&data, &info, sor_steps);

        timings.time[4] += MPI_Wtime();
    }
    timings.time[0] /= max_steps;
    timings.time[4] /= max_steps;
    if(info.rank == 0)
    {
        printf("On p=%d, N=%d, P=%d, Q=%d, nju=%d\n", info.size, info.N, info.P, info.Q, sor_steps);
        printf("Time needed for preconditioned CG: %2.10lf\n", timings.time[4] - timings.time[0]);
    }

    timings.time[0] = 0;
    timings.time[4] = 0;
    for(steps = 0; steps < max_steps; ++steps)
    {
        timings.time[0] += MPI_Wtime();

        CG_solver(&data, &info);

        timings.time[4] += MPI_Wtime();
    }
    timings.time[0] /= max_steps;
    timings.time[4] /= max_steps;
    if(info.rank == 0)
    {
        printf("Time needed for CG               : %2.10lf\n", timings.time[4] - timings.time[0]);
    }

    if(info.rank == 0)
    {
        //printf("\nEND:\n");
        //print_Matrix1D(&data.x);
        //printf("\n");
        //qq(&data.x, info.h);
        //printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    struct Timing t;
//
//    ran(&data.x);
    //zero_up(&data.x);
    //zero_up(&data.r);
    //initialize_starting_values(&data.x, &info);
    //solve(&data.x, &data.r, &info, &t, 1000000);
////////
    //if(info.rank == 0)
    //    print_Matrix1D(&data.x);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info.rank == 1)
//        print_Matrix1D(&data.x);


//
//    printf("\n");
//
//    MPI_Barrier(MPI_COMM_WORLD);
    //print_Matrix1D(&data.s);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(info.rank == 1)
//    {
//        printf("\n");
//        print_Matrix1D(&data.x);
//
//    }
    //solve(&u, &info, &timings);

//    MPI_Barrier(MPI_COMM_WORLD);
//    timings.time[4] = MPI_Wtime();
//
//    if (info.rank == 0)
//    {
//        //print_Matrix1D(&u);
//        printf("\n");
//        //qq(&u, info.h);
//        printf("running on N=%d, P=%d, Q=%d,  p=%d\n", info.N, info.P, info.Q, info.size);
//        printf("start up time %4.5lf, \n", timings.time[1] - timings.time[0]);
//        printf("solve time    %4.5lf, \n", timings.time[4] - timings.time[1]);
//        printf("norm time     %4.5lf, \n", timings.time[3] - timings.time[2]);
//    }


    // Test for valid resutl.

    bool check = 1;//testCompare( &data.x, &info, info.h);
    bool allChecks = 0;
  //  printf("%d %d\n", info.rank, (int)check);
    int sum_test = testCompare(&data.x, &info, info.h);
    int sum_all_test = 0;

    MPI_Reduce(&sum_test, &sum_all_test, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&check, &allChecks, 1, MPI_LOGICAL, MPI_LAND, 1, MPI_COMM_WORLD);
    if(info.rank == 0)
    {
        printf("\nTest results of N = %d and processors = %d: %d\n\n", info.N, info.size, sum_all_test == info.size);
    }

    //free_Matrix1D(&u);

    //printf("=-----------\n%d\n------------", info.rank);
    free_Matrix1D(&data.x);
    free_Matrix1D(&data.s);
    free_Matrix1D(&data.r);
    free_Matrix1D(&data.As);
}

