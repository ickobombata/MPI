#include <mpi.h>
#include <stdio.h>

int main (int argc, char* argv[])
{
    int rank, size;

    MPI_Init (&argc, &argv);      /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */

    printf("Program running on %d processes, I am no process. %d\n", size, rank);

    MPI_Finalize();
    return 0;
}
