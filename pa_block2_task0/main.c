//#include <mpi.h>
#include <stdio.h>
#include "/usr/include/mpi/mpi.h"
#include "SSOR.h"
//#include "matrix_allocation.h"



int main (int argc, char* argv[])
{
    int error;
    struct Info info = {24, 2, 3, 0, 1.0};   //N=24, 48, 120, P = 2, Q = 3, f ≡ 0, ω = 1,
    info.h = 1.0 / (double)(info.N + 1);
    info.l = info.N / info.P;
    info.m = info.N / info.Q;

    MPI_Init (&argc, &argv);      /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &info.rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &info.size);        /* get number of processes */

    //struct Matrix1D u;
    //if(!allocate_Matrix1D(info.l+2, info.m+2, &u))
    //    return; // ERROR here
    //initialize_starting_values(&u);

    //set_up_grid(&info);



    int rank, size, i;
    MPI_Datatype type, type2;
    int q = 30;
    int buffer[30];
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Type_contiguous(1, MPI_INT, &type2);
    MPI_Type_commit(&type2);
    MPI_Type_vector(3, 6, 9, type2, &type);
    MPI_Type_commit(&type);

    if (rank == 0)
    {
        for (i=0; i<q; i++)
            buffer[i] = i;
        MPI_Send(buffer, 1, type, 1, 123, MPI_COMM_WORLD);
    }

    if (rank == 1)
    {
        for (i=0; i<q; i++)
            buffer[i] = -1;
        MPI_Recv(buffer, 1, type, 0, 123, MPI_COMM_WORLD, &status);
        for (i=0; i<q; i++)
            printf("buffer[%d] = %d\n", i, buffer[i]);
        fflush(stdout);
    }

//
//    if(info.rank == 2)
//    {
//        printf("l =%d, m=%d \n", info.l, info.m);
//        int i, j;
//        for(i = 0 ;i < u.n; ++i)
//        {
//            for(j = 0 ; j < u.m; ++j)
//                printf("%2.2lf ,", (*get_Matrix1D(&u, i, j)));
//            printf("\n");
//        }
//
//        printf("%d, %d\n", info.coords[0], info.coords[1]);
//    }

    /// need to pass b also but for now is 0
    //solve(&u, &info);

    //free_Matrix1D(&u);

    MPI_Finalize();
    return 0;
}
