#include <stdio.h>
#include <stdlib.h>
#include "CG.h"
//#include "/usr/include/mpi/mpi.h"

int main(int argc, char* argv[])
{
    MPI_Init (&argc, &argv);

    cg_settup();

    MPI_Finalize();
    return 0;
}
