#include "/usr/include/mpi/mpi.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

void get_user_input(int* arr, size_t size, int rank);
void calculate_subsum(int rank, int* attr, double* sum, int mandatory, int leftover_from, int leftover);
double func(double x);
void fan_in(double* local, int rank, int n, MPI_Status* status);

int main (int argc, char* argv[])
{
    int rank, size, ierr, i;
    double t1, t2;
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    int attr[3];// attr[0] = a , attr[1] = b, attr[2] = n
    double local = 0;

    t1 = MPI_Wtime();
    get_user_input(attr, 3, rank);

    int for_everyone = (int)(((attr[2] << 1) + 1) / size);
    int for_few = ((attr[2] << 1) + 1) - size * for_everyone;
    int taken = size*for_everyone;

    calculate_subsum(rank, attr, &local, for_everyone, taken, for_few);
    //ierr = MPI_Recv(attr, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    //printf("Process %d got %d, %d and  n = %d calculated %lf \n", rank, attr[0], attr[1], attr[2], local);

    fan_in(&local, rank, size, &status);

    if(rank == 0)
    {
        printf("Simpsons Rule for [%d, %d ] with %d subintervals\n", attr[0], attr[1], attr[2]);
        printf("\tyielded the approximation %lf\n", local);
        printf("and took %lf seconds\n", MPI_Wtime() - t1); // measure time
    }

    MPI_Finalize();
    return 0;
}

void fan_in(double* local, int rank, int size, MPI_Status* status)
{
    int N = log2(size);
    int bound = N;
    int k, i;
    int nice_k = 0;
    double received = 1;
    for(k = 1;k <= N; k++)
    {
        if(rank == 0)
        {
            //printf("recv now with rank = %d and k = %d\n", rank, k);
            MPI_Recv(&received, 1, MPI_DOUBLE, (1 << (k-1)), 1, MPI_COMM_WORLD, status);
                *local = *local + received;
        }
        else
        {
            if(rank % (1 << (k-1)) == 0 && rank % (1 << k) != 0) // me == 2i+1
            {
                //printf("send now %lf with rank = %d and k = %d\n",*local, rank, k);
                MPI_Send(local, 1, MPI_DOUBLE, rank - (1 << (k-1)), 1, MPI_COMM_WORLD);
            }
            if(rank % (1 << k) == 0)
            {
                //printf("recv now with rank = %d and k = %d\n", rank, k);
                MPI_Recv(&received, 1, MPI_DOUBLE, rank + (1 << (k-1)), 1, MPI_COMM_WORLD, status);
                *local = *local + received;
            }
        }
    }
}

void calculate_subsum(int rank, int* attr, double* sum, int mandatory, int leftover_from, int leftover)
{
    //printf("processor %d, mandatory = %d \n", rank, mandatory);
    int bound = mandatory * rank + mandatory;
    int i;
    double h = (double)(attr[1] - attr[0]) / (double)attr[2]; // ---> (b-a)/n
    double x_i;
    int mulyiplier = 2;
    for(i = mandatory * rank;i < bound;i++)
    {
        x_i = attr[0] + (i / 2.0) * h; // ----> a + i*h
        if(i == 0 || i == attr[2]<<1) // get rid of this
            mulyiplier = 1;
        else if(i & 1)
            mulyiplier = 4;
        else
            mulyiplier = 2;
        *sum += mulyiplier * func(x_i);
        //printf("processor %d, x_i = %lf mul = %d func = %lf\n", rank, x_i, mulyiplier, func(x_i));
    }

    if(rank < leftover)
    {
        //printf("processor %d, leftover = %d \n", rank, leftover_from);

        x_i = attr[0] + ((leftover_from + rank)/2.0) * h;
        if(rank == 0 || rank == attr[2]<<1) // get rid of this
            mulyiplier = 1;
        else if(rank & 1)
            mulyiplier = 4;
        else
            mulyiplier = 2;
        *sum += mulyiplier * func(x_i);
    }
}

double func(double x)
{
	// double x2 = x*x
	// double x4 = x2*x2
	// double result = x4*x2;
    return pow(x, 6);
}

void get_user_input(int* arr, size_t size,int rank)
{
    int ierr;
    if(rank == 0)
    {
        printf("Enter the interval bounds (a, b): \n");
        ierr = scanf("%d", &arr[0]);
        ierr = scanf("%d", &arr[1]);
        printf("Enter the number of subintervals(n): \n");
        ierr = scanf("%d", &arr[2]);
    }
    ierr = MPI_Bcast(arr, 3, MPI_INT, 0, MPI_COMM_WORLD);
}
