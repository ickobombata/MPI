#include "/usr/include/mpi/mpi.h"
#include <stdio.h>
#include <string.h>

int main (int argc, char* argv[])
{
    int rank, size, ierr, i;
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);


    double number;
    if(rank == 0)
    {
        printf("Enter number: ");
        ierr = scanf("%lf", &number);
        printf("Process 0 got %lf and computed %lf\n", number, number + number);
        if(size > 1)
        {
            number = number + number;
            ierr = MPI_Send(&number, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
    }
    else if (rank == size - 1)
    {
        ierr = MPI_Recv(&number, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        printf("Process %d got %lf and computed %lf\n", rank, number, number + number);
    }
    else
    {
        ierr = MPI_Recv(&number, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        printf("Process %d got %lf and computed %lf\n", rank, number, number + number);
        number = number + number;
        ierr = MPI_Send(&number, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
